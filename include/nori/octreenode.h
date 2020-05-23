#pragma once

#include <nori/object.h>
#include <nori/common.h>
#include <nori/mesh.h>
NORI_NAMESPACE_BEGIN
class OcTreeNode {
private:
    BoundingBox3f boundingBox;
    bool isLeaf;
    std::vector<OcTreeNode*> children;
    const Mesh* pMesh;
    std::vector<int> triangleIndices;
public:
    OcTreeNode(const BoundingBox3f& bb);
    OcTreeNode();
    ~OcTreeNode();
    
    void SetBoundingBox(const BoundingBox3f& bb);
    const BoundingBox3f& GetBoundingBox() const;
    bool rayIntersect(const Ray3f &ray, Intersection &its, bool shadowRay) const;
private:
    bool ShouldSubdivide(int maxDepth) const;
    void SubdivideRecursive(int maxDepth);
    void SetMeshTriangles(const Mesh* mesh, const std::vector<int>& triangles);
    friend OcTreeNode* BuildOctree(const Mesh* mesh, int maxDepth);
};

OcTreeNode* BuildOctree(const Mesh* mesh, int maxDepth);

NORI_NAMESPACE_END
