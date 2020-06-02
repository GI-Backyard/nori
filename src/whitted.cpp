#include <nori/integrator.h>
#include <nori/scene.h>
#include <algorithm>

NORI_NAMESPACE_BEGIN

class WhittedIntegrator : public Integrator {
public:
    WhittedIntegrator(const PropertyList& props) {
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        Intersection its;
        if(!scene->rayIntersect(ray, its)) {
            return Color3f(0.0, 0.0, 0.0);
        }
        Normal3f n = its.shFrame.n.cwiseAbs();
        return Color3f(n.x(), n.y(), n.z());
    }
    /// Return a human-readable description for debugging purposes
    std::string toString() const {
        return "PointLightingIntegrator[]";
    }
};
NORI_REGISTER_CLASS(WhittedIntegrator, "whitted")
NORI_NAMESPACE_END
