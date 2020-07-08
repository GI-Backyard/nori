#include <nori/integrator.h>
#include <nori/scene.h>
#include <algorithm>
#include <nori/emitter.h>
#include <nori/warp.h>
#include <nori/bsdf.h>

NORI_NAMESPACE_BEGIN

class WhittedIntegrator : public Integrator {
public:
    WhittedIntegrator(const PropertyList& props) {
    }

    Color3f EvaluateDiffuse(const Scene *scene, Sampler *sampler, const Intersection& its, const Ray3f &lastRay) const
    {
        Color3f Ltotal(0, 0, 0);
        // Le
        if(its.mesh->getEmitter())
        {
            EmitterQueryRecord record;
            Ltotal = its.mesh->getEmitter()->eval(record);
        }
        
        // sampling Light
        const std::vector<Emitter*>& lights = scene->getLights();
        uint32_t sampledLight = static_cast<uint32_t>(sampler->next1D() * lights.size());
        EmitterQueryRecord lightRecord;
        Color3f Le = lights[sampledLight]->sample(lightRecord, Point3f(sampler->next1D(), sampler->next1D(), sampler->next1D()));
        float pdf = lights[sampledLight]->pdf(lightRecord);
        pdf = pdf / lights.size();
        
        Vector3f wo = -lastRay.d;
        wo.normalize();
        Vector3f wi = lightRecord.position - its.p;
        float lengthwi = wi.norm();
        wi.normalize();
        
        Ray3f secondaRay(its.p, wi);
        secondaRay.maxt = lengthwi - Epsilon;
        if(!scene->rayIntersect(secondaRay))
        {
			float nDotL = its.shFrame.n.dot(wi);
			float lightNDotNegL = lightRecord.normal.dot(-wi);
			float geom = its.shFrame.n.dot(wi) * lightRecord.normal.dot(-wi) / (lengthwi * lengthwi);
			if (nDotL < 0 || lightNDotNegL < 0)
			{
				geom = 0;
			}
			else
			{
				// bsdf
				Color3f reflection = its.mesh->getBSDF()->eval(BSDFQueryRecord(its.toLocal(wi), its.toLocal(wo), EMeasure::ESolidAngle));
				Ltotal += Le * reflection * geom / pdf;
			}

        }
        
        
        return Ltotal;
    }
    
    // note: this could be iterative for multiple specular
    Color3f EvaluateSpecular(const Scene *scene, Sampler *sampler, const Intersection& its, const Ray3f &lastRay) const
    {
        const BSDF* bsdf = its.mesh->getBSDF();
        BSDFQueryRecord record(its.toLocal(-lastRay.d));
        Color3f bsdfdivPdf = bsdf->sample(record, sampler->next2D());
        Ray3f secondaRay(its.p, its.toWorld(record.wo));
        Intersection nextIts;
        const float earlyExitThreshold = 0.95f;
        if(scene->rayIntersect(secondaRay, nextIts))
        {
            if(sampler->next1D() < earlyExitThreshold)
            {
                return 1.0/0.95 * bsdfdivPdf * EvaluateWhitted(scene, sampler, nextIts, secondaRay);
            }
            else
            {
                return Color3f(0, 0, 0);
            }
        }
        else
        {
            return Color3f(0, 0, 0);
        }
    }
    
    Color3f EvaluateWhitted(const Scene *scene, Sampler *sampler, const Intersection& its, const Ray3f &lastRay) const
    {
        if(!its.mesh->getBSDF()->isDiffuse())
        {
            return EvaluateSpecular(scene, sampler, its, lastRay);
        }
        return EvaluateDiffuse(scene, sampler, its, lastRay);
    }
    
    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        Intersection its;
        if(!scene->rayIntersect(ray, its)) {
            return Color3f(0.1f, 0.1f, 0.1f);
        }
        return EvaluateWhitted(scene, sampler, its, ray);
        
    }
    /// Return a human-readable description for debugging purposes
    std::string toString() const {
        return "PointLightingIntegrator[]";
    }
};
NORI_REGISTER_CLASS(WhittedIntegrator, "whitted")
NORI_NAMESPACE_END
