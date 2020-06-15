#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>
#include <nori/bsdf.h>
NORI_NAMESPACE_BEGIN

void CheckColor(const std::string& name, Color3f& color)
{
	if (color.x() <= 0 || color.y() <= 0 || color.z() <= 0)
	{
		std::cout << name << " (" << color.x() << "," << color.y() << "," << color.z() << ")" << std::endl;
	}
}

float maxComponent(Color3f& color)
{
	float result = color.x();
	if (color.y() > result) result = color.y();
	if (color.z() > result) result = color.z();
	return result;
}

#define CHECK_COLOR(a) CheckColor(#a, a);
class PathMatsIntegrator : public Integrator {
public:
    PathMatsIntegrator(const PropertyList& props) {
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        Color3f Ltotal(0, 0, 0);
        Intersection its;
		Ray3f tracedRay = ray;
		Color3f rayThrough = Color3f(1.0f);
        const int MAX_BOUNCE = 8;
		float accumunatedEta = 1.0f;
        for(int bounce = 0; bounce < MAX_BOUNCE; ++bounce)
        {
			bool intersected = scene->rayIntersect(tracedRay, its);
			{
				if (intersected && its.mesh->getEmitter())
				{
					EmitterQueryRecord record;
					Ltotal += rayThrough * its.mesh->getEmitter()->eval(record);
				}
			}

			if (!intersected)
				break;

			Vector3f wo = -tracedRay.d;
			wo.normalize();

			// sample bsdf
			const BSDF* bsdf = its.mesh->getBSDF();
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
        return "PathMatsIntegrator[]";
    }
};
NORI_REGISTER_CLASS(PathMatsIntegrator, "path_mats")
NORI_NAMESPACE_END
