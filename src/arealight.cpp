
#include <nori/emitter.h>
#include <nori/mesh.h>

NORI_NAMESPACE_BEGIN

class AreaLight : public Emitter
{
public:
    virtual Color3f sample(EmitterQueryRecord& record, const Point3f& sample) const
    {
        if(geometry)
        {
            Point3f position;
            Vector3f normal;
            float pdf;
            geometry->uniformSample(sample, position, normal, pdf);
            record.position = position;
            record.normal = normal;
        }
        return radiance;
    }

    virtual Color3f eval(const EmitterQueryRecord& record) const
    {
        return radiance;
    }

    virtual float pdf(const EmitterQueryRecord& record) const
    {
        return geometry != nullptr ? geometry->getFacePdf()->getNormalization() : 1.0f;
    }

    virtual std::string toString() const
    {
        return "AreaLight[]";
    }

    AreaLight(const PropertyList& props)
        : geometry(nullptr)
    {
        radiance = props.getColor("radiance");
    }

    void setParent(NoriObject *parent)
    {
        Emitter::setParent(parent);
        if(parent->getClassType() == EClassType::EMesh)
        {
            geometry = dynamic_cast<const Mesh*>(parent);
        }
    }
private:
    Color3f radiance;
    const Mesh* geometry;
};

NORI_REGISTER_CLASS(AreaLight, "area")

NORI_NAMESPACE_END
