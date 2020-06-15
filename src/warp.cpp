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

#include <nori/warp.h>
#include <nori/vector.h>
#include <nori/frame.h>
#include <nori/common.h>
NORI_NAMESPACE_BEGIN

Point2f Warp::squareToUniformSquare(const Point2f &sample) {
    return sample;
}

float Warp::squareToUniformSquarePdf(const Point2f &sample) {
    return ((sample.array() >= 0).all() && (sample.array() <= 1).all()) ? 1.0f : 0.0f;
}

Point2f Warp::squareToTent(const Point2f &sample) {
    auto tent1D = [](float delta)
    {
        if(delta < 0 && delta > 1) return 0.f;
        else if(delta >=0 && delta <= 0.5) return sqrt(2 * delta) - 1.0f;
        else return 1 - sqrt(2.0f - 2 * delta);
    };
    return Point2f(tent1D(sample.x()), tent1D(sample.y()));
}

float Warp::squareToTentPdf(const Point2f &p) {
    // the area is [-1, 1] * [-1, 1] square
    return (abs(p.x()) <= 1.0f && abs(p.y()) <= 1.0f) ? (1 - abs(p.x())) * ( 1 - abs(p.y())) : 0.0f;
}

Point2f Warp::squareToUniformDisk(const Point2f &sample) {
#define NAIIVE_UNIFORMDISK_DISTRIBUTION 0
#if NAIIVE_UNIFORMDISK_DISTRIBUTION
    //    naiive implementation
    float r = sqrt(sample.x());
    float theta = 2 * M_PI * sample.y();
    return Point2f(r * cos(theta), r * sin(theta));
#else
    const float x = sample.x() * 2 - 1.0f;
    const float y = sample.y() * 2 - 1.0f;
    if(x == 0 && y == 0) return Point2f(0, 0);
    
    float r(0), theta(0);
	const float piDivide4 = M_PI * 0.25;
    if(x > 0 && abs(y) <= x)
    {
        r = x;
        theta = y / x * M_PI * 0.25;
    }
    else if( x < 0 && abs(y) <= -x)
    {
        r = -x;
        theta = M_PI - y / x * M_PI * 0.25;
    }
    else if (y > 0 && abs(x) < y)
    {
        r = y;
        theta = M_PI_2 - x / y * M_PI * 0.25;
    }
    else
    {
        r = -y;
        theta = 3 * M_PI_2 + x / y * M_PI * 0.25;
    }
    return Point2f(r * cos(theta), r * sin(theta));
#endif
}

float Warp::squareToUniformDiskPdf(const Point2f &p) {
    return (p.x() * p.x() + p.y() * p.y() <= 1) ? INV_PI : 0.0f;
}

Vector3f Warp::squareToUniformSphere(const Point2f &sample) {
    float costheta = 1 - 2 * sample.x();
    float sintheta = sqrt( 1 - costheta * costheta);
    float phi = 2 * M_PI * sample.y();
    return Vector3f(sintheta * cos(phi), sintheta * sin(phi), costheta);
}

float Warp::squareToUniformSpherePdf(const Vector3f &v) {
    return 0.25 * INV_PI;
}

Vector3f Warp::squareToUniformHemisphere(const Point2f &sample) {
    float costheta = sample.y();
    float sintheta = sqrt( 1 - costheta * costheta);
    float phi = 2 * M_PI * sample.x();
    return Vector3f(sintheta * cos(phi), sintheta * sin(phi), costheta);
}

float Warp::squareToUniformHemispherePdf(const Vector3f &v) {
    return (v.z() >= 0) ? (0.5 * INV_PI) : 0.0f;
}

Vector3f Warp::squareToCosineHemisphere(const Point2f &sample) {
    float costheta = sqrt( 1 - sample.x());
    float sintheta = sqrt(sample.x());
    float phi = 2 * M_PI * sample.y();
    return Vector3f(sintheta * cos(phi), sintheta * sin(phi), costheta);
}

float Warp::squareToCosineHemispherePdf(const Vector3f &v) {
    return (v.z() >= 0) ? (v.z() * INV_PI) : 0.0f;
}

Vector3f Warp::squareToBeckmann(const Point2f &sample, float alpha) {
//    float tanthetaSquared = - alpha * alpha * log(sample.x());
//    float costheta = sqrt(1/(1 + tanthetaSquared));
//    float sintheta = sqrt(tanthetaSquared/(1 + tanthetaSquared));
    float costheta = sqrt(1/ (1 - alpha * alpha * log(sample.x())));
    float sintheta = sqrt( 1 - costheta * costheta);
    float phi = 2 * M_PI * sample.y();
    return Vector3f(sintheta * cos(phi), sintheta * sin(phi), costheta);
}

float Warp::squareToBeckmannPdf(const Vector3f &m, float alpha) {
    float costheta = m.z();
    float sintheta = sqrt( 1 - m.z() * m.z());
    float tantheta = sintheta / costheta;
    if(m.z() > 0) return INV_PI * exp(-tantheta * tantheta / (alpha * alpha)) / (alpha * alpha * pow(costheta, 3));
    else return 0.0f;
}

NORI_NAMESPACE_END
