#include <nori/integrator.h>
#include <nori/scene.h>
#include <algorithm>
#include <nori/warp.h>
#include <pcg32.h>

NORI_NAMESPACE_BEGIN

class AOIntegrator : public Integrator {
public:
    AOIntegrator(const PropertyList& props) {
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        Intersection its;
        if(!scene->rayIntersect(ray, its)) {
            return Color3f(0.2, 0.2, 0.2);
        }
        float aoIntegral = 0;
        Vector3f dir = Warp::squareToCosineHemisphere(sampler->next2D());
        dir = its.toWorld(dir);
        dir.normalize();
        if(!scene->rayIntersect(Ray3f(its.p, dir))) aoIntegral = aoIntegral + 1;
        return Color3f(aoIntegral, aoIntegral, aoIntegral);
    }
    /// Return a human-readable description for debugging purposes
    std::string toString() const {
        return "AOIntegrator[]";
    }
    Point3f position;
    Color3f energy;
};
NORI_REGISTER_CLASS(AOIntegrator, "ao")
NORI_NAMESPACE_END
