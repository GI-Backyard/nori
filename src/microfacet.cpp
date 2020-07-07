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
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class Microfacet : public BSDF {
public:
    Microfacet(const PropertyList &propList) {
        /* RMS surface roughness */
        m_alpha = propList.getFloat("alpha", 0.1f);

        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);

        /* Albedo of the diffuse base material (a.k.a "kd") */
        m_kd = propList.getColor("kd", Color3f(0.5f));

        /* To ensure energy conservation, we must scale the 
           specular component by 1-kd. 

           While that is not a particularly realistic model of what 
           happens in reality, this will greatly simplify the 
           implementation. Please see the course staff if you're 
           interested in implementing a more realistic version 
           of this BRDF. */
        m_ks = 1 - m_kd.maxCoeff();
    }

    /// Evaluate the BRDF for the given pair of directions
    Color3f eval(const BSDFQueryRecord &bRec) const {
		Color3f result = m_kd * INV_PI;
		//Vector3f wi = bRec.wi;
		//wi.x() = abs(wi.x());
		Vector3f wh = bRec.wi + bRec.wo;
		wh.normalize();
		float g = G1(bRec.wi, wh, Vector3f(0, 0, 1)) * G1(bRec.wo, wh, Vector3f(0, 0, 1));
		// todo : verify d calculation
		float d = 1.0f;
		float costheI = clamp(bRec.wi.z(), 0.0f, 1.0f);
		float costheO = clamp(bRec.wo.z(), 0.0f, 1.0f);
		float costheH = clamp(wh.z(), 0.0f, 1.0f);
		float f = fresnel(clamp(bRec.wi.dot(wh), 0.0f, 1.0f), bRec.eta, 1.0f);
		result += m_ks * d * f * g / (4.0f * costheI * costheO * costheH);
		if (!result.isValid())
		{
			cerr << "Integrator: computed an invalid radiance value: " << result.toString() << endl;
		}
		return result;
    }

    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const {
		Vector3f wh = bRec.wi + bRec.wo;
		wh.normalize();
		return Warp::squareToBeckmannPdf(wh, m_alpha);
    }

    /// Sample the BRDF
    Color3f sample(BSDFQueryRecord &bRec, const Point2f &_sample) const {
		// sample wh
		Vector3f wh = Warp::squareToBeckmann(_sample, m_alpha);
		bRec.wo = wh * 2 - bRec.wi;
		bRec.eta = m_extIOR / m_intIOR;
		bRec.measure = EMeasure::ESolidAngle;

        // Note: Once you have implemented the part that computes the scattered
        // direction, the last part of this function should simply return the
        // BRDF value divided by the solid angle density and multiplied by the
        // cosine factor from the reflection equation, i.e.
        return eval(bRec) * Frame::cosTheta(bRec.wo) / pdf(bRec);
    }

    bool isDiffuse() const {
        /* While microfacet BRDFs are not perfectly diffuse, they can be
           handled by sampling techniques for diffuse/non-specular materials,
           hence we return true here */
        return true;
    }

    std::string toString() const {
        return tfm::format(
            "Microfacet[\n"
            "  alpha = %f,\n"
            "  intIOR = %f,\n"
            "  extIOR = %f,\n"
            "  kd = %s,\n"
            "  ks = %f\n"
            "]",
            m_alpha,
            m_intIOR,
            m_extIOR,
            m_kd.toString(),
            m_ks
        );
    }
private:
	float G1(const Vector3f& wv, const Vector3f& wh, const Vector3f& n) const
	{
		float lamda = wv.dot(wh) / wv.dot(n);
		lamda = lamda <= 0 ? 0.f : 1.0f;
		float costheV = wv.dot(n);
		float tantheV = sqrt(1 - costheV * costheV) / costheV;
		float b = 1 / (m_alpha * tantheV);
		b = b < 1.6f ? (3.535f * b + 2.181f * b * b) / (1.f + 2.276 * b + 2.577 * b * b) : 1.0f;
		return lamda * b;
	}
private:
    float m_alpha;
    float m_intIOR, m_extIOR;
    float m_ks;
    Color3f m_kd;
};

NORI_REGISTER_CLASS(Microfacet, "microfacet");
NORI_NAMESPACE_END
