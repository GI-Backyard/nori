#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>
#include <nori/bsdf.h>
NORI_NAMESPACE_BEGIN

static void CheckColor(const std::string& name, Color3f& color)
{
	if (color.x() <= 0 || color.y() <= 0 || color.z() <= 0)
	{
		std::cout << name << " (" << color.x() << "," << color.y() << "," << color.z() << ")" << std::endl;
	}
}

static float maxComponent(Color3f& color)
{
	float result = color.x();
	if (color.y() > result) result = color.y();
	if (color.z() > result) result = color.z();
	return result;
}

inline float BalanceHeuristic(int nf, float fPdf, int ng, float gPdf) {
    return (nf * fPdf) / (nf * fPdf + ng * gPdf);
}

inline float PowerHeuristic(int nf, float fPdf, int ng, float gPdf) {
    float f = nf * fPdf, g = ng * gPdf;
    return (f * f) / (f * f + g * g);
}

Color3f EMSSampleLights(const Scene *scene, const Ray3f& ray, const Intersection& its, Sampler* sampler)
{
	Color3f result(0.f);
	const std::vector<Emitter*>& lights = scene->getLights();
	uint32_t sampledLight = sampler->next1D() * lights.size();
	// mis sample light
	{
		EmitterQueryRecord lightRecord;
		Color3f Le = lights[sampledLight]->sample(lightRecord, Point3f(sampler->next1D(), sampler->next1D(), sampler->next1D()));
		float lightPdf = lights[sampledLight]->pdf(lightRecord);

		Vector3f wo = -ray.d;
		wo.normalize();
		Vector3f wi = lightRecord.position - its.p;
		float lengthwi = wi.norm();
		wi.normalize();

		Ray3f shadowRay(its.p, wi);
		shadowRay.maxt = lengthwi - Epsilon;
		if (!scene->rayIntersect(shadowRay))
		{
			float geom = its.shFrame.n.dot(wi) * lightRecord.normal.dot(-wi) / (lengthwi * lengthwi);
			if (geom < 0) geom = 0;
			// bsdf
			Color3f reflection = its.mesh->getBSDF()->eval(BSDFQueryRecord(its.toLocal(wo), its.toLocal(wi), EMeasure::ESolidAngle));
			float bsdfPdf = its.mesh->getBSDF()->pdf(BSDFQueryRecord(its.toLocal(wo), its.toLocal(wi), EMeasure::ESolidAngle));
			result = Le * reflection * geom / lightPdf;

		}
	}

	return result;
}

Color3f MISSampleLights(const Scene *scene, const Ray3f& ray, const Intersection& its, Sampler* sampler)
{
	Color3f result(0.f);
	const std::vector<Emitter*>& lights = scene->getLights();
	uint32_t sampledLight = sampler->next1D() * lights.size();
	// mis sample light
	{
		EmitterQueryRecord lightRecord;
		Color3f Le = lights[sampledLight]->sample(lightRecord, Point3f(sampler->next1D(), sampler->next1D(), sampler->next1D()));
		float lightPdf = lights[sampledLight]->pdf(lightRecord);

		Vector3f wo = -ray.d;
		wo.normalize();
		Vector3f wi = lightRecord.position - its.p;
		float lengthwi = wi.norm();
		wi.normalize();

		Ray3f shadowRay(its.p, wi);
		shadowRay.maxt = lengthwi - Epsilon;
		if (!scene->rayIntersect(shadowRay))
		{
			float geom = its.shFrame.n.dot(wi) * lightRecord.normal.dot(-wi) / (lengthwi * lengthwi);
			if (geom < 0) geom = 0;
			// bsdf
			Color3f reflection = its.mesh->getBSDF()->eval(BSDFQueryRecord(its.toLocal(wo), its.toLocal(wi), EMeasure::ESolidAngle));
			float bsdfPdf = its.mesh->getBSDF()->pdf(BSDFQueryRecord(its.toLocal(wo), its.toLocal(wi), EMeasure::ESolidAngle));
			float weight = PowerHeuristic(1, lightPdf, 1, bsdfPdf);
			result = Le * reflection * geom * weight / lightPdf;
			
		}
	}

	// mis sample bsdf
	{
		Vector3f wo = -ray.d;
		wo.normalize();
		BSDFQueryRecord bsdfRecord(its.toLocal(wo));
		Color3f fdivPdf = its.mesh->getBSDF()->sample(bsdfRecord, sampler->next2D());
		float bsdfPdf = its.mesh->getBSDF()->pdf(bsdfRecord);
		float weight = 1.0f;
		float lightPdf = 0.0f;
		Color3f Le(0.f);
		Ray3f shadowRay(its.p, its.toWorld(bsdfRecord.wo));
		Intersection shadowIntersection;
		bool intersected = scene->rayIntersect(shadowRay, shadowIntersection);
		if (intersected && shadowIntersection.mesh->getEmitter() == lights[sampledLight] && bsdfPdf != 0)
		{
			EmitterQueryRecord lightRecord;
			lightRecord.position = shadowIntersection.p;
			lightRecord.normal = shadowIntersection.geoFrame.n;
			lightPdf = lights[sampledLight]->pdf(lightRecord);
			weight = PowerHeuristic(1, bsdfPdf, 1, lightPdf);
			Le = lights[sampledLight]->eval(lightRecord);
		}

		result += Le * fdivPdf * weight;
	}

	return result * lights.size();
}

class PathMISIntegrator : public Integrator {
public:
    PathMISIntegrator(const PropertyList& props) {
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        Color3f Ltotal(0, 0, 0);
        Intersection its;
		Ray3f tracedRay = ray;
		Color3f rayThrough = Color3f(1.0f);
        const int MAX_BOUNCE = 10;
		float accumunatedEta = 1.0f;
		bool isSpecular = false;
        for(int bounce = 0; bounce < MAX_BOUNCE; ++bounce)
        {
			bool intersected = scene->rayIntersect(tracedRay, its);
			{
				// for specular and eye ray, emitter need to be calculated seperately
				if (intersected && its.mesh->getEmitter() && (bounce == 0 || isSpecular))
				{
					EmitterQueryRecord record;
					Ltotal += rayThrough * its.mesh->getEmitter()->eval(record);
				}
			}

			if (!intersected)
				break;

			// calculate lighting
            Ltotal += rayThrough * MISSampleLights(scene, tracedRay, its, sampler);

			Vector3f wo = -tracedRay.d;
			wo.normalize();

			// sample bsdf
			const BSDF* bsdf = its.mesh->getBSDF();
			isSpecular = !bsdf->isDiffuse();
			BSDFQueryRecord record(its.toLocal(wo));

			Color3f f = bsdf->sample(record, sampler->next2D());
			if (f.nonZeros() == 0) break;
			//  next ray
			tracedRay = Ray3f (its.p, its.toWorld(record.wo));

			rayThrough = rayThrough * f;
			
			// russian roulette
			accumunatedEta *= record.eta;
			if (bounce > 3)
			{
				float continueProbility = maxComponent(rayThrough) * accumunatedEta * accumunatedEta;
				if (continueProbility > 0.99f) continueProbility = 0.99f;
				if (sampler->next1D() > continueProbility) break;
				else rayThrough /= continueProbility;

			}			

        }

		return Ltotal;
        //if(!scene->rayIntersect(ray, its)) {
        //    return Color3f(0, 0, 0);
        //}
        //Normal3f n = its.shFrame.n.cwiseAbs();
        //return Color3f(n.x(), n.y(), n.z());
    }
    /// Return a human-readable description for debugging purposes
    std::string toString() const {
        return "PathMISIntegrator[]";
    }
};
NORI_REGISTER_CLASS(PathMISIntegrator, "path_mis")
NORI_NAMESPACE_END
