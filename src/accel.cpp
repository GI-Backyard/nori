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

#include <nori/accel.h>
#include <Eigen/Geometry>

NORI_NAMESPACE_BEGIN

void Accel::addMesh(Mesh *mesh) {
    // if (m_mesh)
    //     throw NoriException("Accel: only a single mesh is supported!");
    if(std::find(m_meshs.begin(), m_meshs.end(), mesh) == m_meshs.end())
    {
        m_meshs.push_back(mesh);
        if(m_meshs.size() == 1) m_bbox = m_meshs[0]->getBoundingBox();
        else m_bbox.expandBy(m_meshs[m_meshs.size() - 1]->getBoundingBox());
    }
    
}

void Accel::build() {
    // build Triangle meshes
    std::vector<Triangle> triangles;
    for(const Mesh* mesh : m_meshs)
    {
        uint32_t currentSize = triangles.size();
        triangles.resize(currentSize + mesh->getTriangleCount());
        for(uint32_t index = 0; index < mesh->getTriangleCount(); ++index)
        {
            triangles[currentSize + index].mesh = mesh;
            triangles[currentSize + index].triangleIndex = index;
        }
    }
    /* build oc tree of max depth 6 */
    m_ocTree = BuildOctree(triangles, m_bbox, 6);
}

bool Accel::rayIntersect(const Ray3f &ray_, Intersection &its, bool shadowRay) const {
    its.t = MAXFLOAT;
    return m_ocTree->rayIntersect(ray_, its, shadowRay);

}

NORI_NAMESPACE_END

