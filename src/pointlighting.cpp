#include <nori/integrator.h>
#include <nori/scene.h>
#include <algorithm>

NORI_NAMESPACE_BEGIN

class PointLightingIntegrator : public Integrator {
public:
    PointLightingIntegrator(const PropertyList& props) {
        position = props.getPoint("position");
        energy = props.getColor("energy");
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        Intersection its;
        if(!scene->rayIntersect(ray, its)) {
            return Color3f(0.0, 0.0, 0.0);
        }
        Vector3f li = position - its.p;
        float length = li.norm();
        li.normalize();
        Ray3f shadowRay(its.p, li);
        float visible = 1.0f;
        if(scene->rayIntersect(shadowRay)) visible = 0;
        float ndotl = std::max(0.0f, li.dot(its.shFrame.n));
        return energy * ndotl * visible / (M_PI * M_PI * 4 * length * length);
//        Normal3f n = its.shFrame.n.cwiseAbs();
//        return Color3f(n.x(), n.y(), n.z());
    }
    /// Return a human-readable description for debugging purposes
    std::string toString() const {
        return "PointLightingIntegrator[]";
    }
    Point3f position;
    Color3f energy;
};
NORI_REGISTER_CLASS(PointLightingIntegrator, "simple")
NORI_NAMESPACE_END
