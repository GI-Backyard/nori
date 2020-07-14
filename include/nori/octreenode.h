#pragma once

#include <nori/object.h>
#include <nori/common.h>
#include <nori/mesh.h>
NORI_NAMESPACE_BEGIN

struct Triangle
{
    const Mesh* mesh;
    int32_t triangleIndex;
};

class OcTreeNode {
private:
    BoundingBox3f boundingBox;
    bool isLeaf;
    std::vector<OcTreeNode*> children;
    std::vector<Triangle> triangles;
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
    void SetMeshTriangles(const std::vector<Triangle>& triangles);
    friend OcTreeNode* BuildOctree(const std::vector<Triangle>& triangles, const BoundingBox3f& bbox, int maxDepth);
};

OcTreeNode* BuildOctree(const std::vector<Triangle>& triangles, const BoundingBox3f& bbox, int maxDepth);

NORI_NAMESPACE_END
