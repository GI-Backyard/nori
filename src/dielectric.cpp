/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/common.h>
NORI_NAMESPACE_BEGIN

/// Ideal dielectric BSDF
class Dielectric : public BSDF {
public:
    Dielectric(const PropertyList &propList) {
        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);
    }

    Color3f eval(const BSDFQueryRecord &) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return Color3f(0.0f);
    }

    float pdf(const BSDFQueryRecord &) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return 0.0f;
    }

    void reflect(BSDFQueryRecord &bRec) const
    {
        // mirror reflect
        // Reflection in local coordinates
        bRec.wo = Vector3f(
                           -bRec.wi.x(),
                           -bRec.wi.y(),
                           bRec.wi.z()
                           );
        bRec.measure = EDiscrete;
        
        /* Relative index of refraction: no change */
        bRec.eta = 1.0f;
    }
    
    bool refract(const Vector3f &wi, const Vector3f& normal,float extI, float extT, BSDFQueryRecord &bRec) const
    {
        float costhetaI = wi.dot(normal);
        float sin2thetaI = std::max(0.0f, 1.0f - costhetaI * costhetaI);
        float sin2thetaT = extI * extI * sin2thetaI / (extT * extT);
        if(sin2thetaT > 1) return false;
        float costhetaT = std::sqrt(1 - sin2thetaT);
        
        bRec.wo = -extI / extT * wi + (extI / extT * costhetaI - costhetaT) * normal;
        bRec.eta = extI / extT;
        bRec.measure = EDiscrete;
        return true;
    }
    
    Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const {
        float extI = m_extIOR;
        float extT = m_intIOR;
        Vector3f normal(0, 0, 1);
        if (Frame::cosTheta(bRec.wi) <= 0)
        {
            normal.z() = -1.f;
            std::swap(extI, extT);
        }
        bool refractSuccess = refract(bRec.wi, normal, extI, extT, bRec);
        if(!refractSuccess)
        {
            reflect(bRec);
            return Color3f(1.0f);
        }
        else
        {
            float frFactor = fresnel(std::abs(bRec.wi.z()), extI, extT);
            if(sample.x() < frFactor)
            {
                reflect(bRec);
                return Color3f(1.0f);
            }
            else
            {
                float r = extT * extT / (extI * extI) * (1.0f - frFactor);
                return Color3f(r);
            }
        }
        
    }

    std::string toString() const {
        return tfm::format(
            "Dielectric[\n"
            "  intIOR = %f,\n"
            "  extIOR = %f\n"
            "]",
            m_intIOR, m_extIOR);
    }
private:
    float m_intIOR, m_extIOR;
};

NORI_REGISTER_CLASS(Dielectric, "dielectric");
NORI_NAMESPACE_END
