#include <nori/integrator.h>
#include <nori/scene.h>

NORI_NAMESPACE_BEGIN

class PathMatsIntegrator : public Integrator {
public:
    PathMatsIntegrator(const PropertyList& props) {
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        Color3f Ltotal(0, 0, 0);
        Intersection its;
        const int MAX_BOUNCE = 8;
        int bounce = 0;
        while (bounce < MAX_BOUNCE)
        {
            /* code */
        }
        
        if(!scene->rayIntersect(ray, its)) {
            return Color3f(0, 0, 0);
        }
        Normal3f n = its.shFrame.n.cwiseAbs();
        return Color3f(n.x(), n.y(), n.z());
    }
    /// Return a human-readable description for debugging purposes
    std::string toString() const {
        return "PathMatsIntegrator[]";
    }
};
NORI_REGISTER_CLASS(PathMatsIntegrator, "path_mats")
NORI_NAMESPACE_END
