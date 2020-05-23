#include <nori/octreenode.h>
#include <Eigen/Geometry>
NORI_NAMESPACE_BEGIN

OcTreeNode* BuildOctree(const Mesh* mesh, int maxDepth)
{
    if(mesh == nullptr || mesh->getTriangleCount() == 0) return nullptr;
    OcTreeNode* root = new OcTreeNode(mesh->getBoundingBox());
    std::vector<int> triangles(mesh->getTriangleCount());
    for(int index = 0; index < mesh->getTriangleCount(); ++index)
    {
        triangles[index] = index;
    }
    root->SetMeshTriangles(mesh, triangles);
    root->SubdivideRecursive(maxDepth);
    return root;
}

OcTreeNode::OcTreeNode(const BoundingBox3f& bb)
: boundingBox(bb)
, isLeaf(true)
{
    
}


OcTreeNode::OcTreeNode()
: isLeaf(true)
{
    
}

OcTreeNode::~OcTreeNode()
{
    if(!isLeaf)
    {
        for(int index = 0; index < children.size(); ++index)
        {
            delete children[index];
        }
        children.clear();
    }
}

void OcTreeNode::SetBoundingBox(const BoundingBox3f& bb)
{
    boundingBox = bb;
}

const BoundingBox3f& OcTreeNode::GetBoundingBox() const
{
    return boundingBox;
}

bool OcTreeNode::ShouldSubdivide(int maxDepth) const
{
    return maxDepth > 0 && triangleIndices.size() > 10;
}

void OcTreeNode::SubdivideRecursive(int maxDepth)
{
    if(ShouldSubdivide(maxDepth))
    {
        isLeaf = false;
        BoundingBox3f bbs[8];
        std::vector<int> subdivideTriangles[8];
        for(int index = 0; index < 8; ++index)
        {
            bbs[index].reset();
            bbs[index].expandBy(boundingBox.getCenter());
            bbs[index].expandBy(boundingBox.getCorner(index));
            for(int tIndex= 0; tIndex < triangleIndices.size(); ++tIndex)
            {
                if(bbs[index].overlaps(pMesh->getBoundingBox(triangleIndices[tIndex])))
                {
                    subdivideTriangles[index].push_back(triangleIndices[tIndex]);
                }
            }
            
            OcTreeNode* child = new OcTreeNode(bbs[index]);
            child->SetMeshTriangles(pMesh, subdivideTriangles[index]);
            child->SubdivideRecursive(maxDepth - 1);
            children.push_back(child);
        }
        
        triangleIndices.clear();
    }
}

bool OcTreeNode::rayIntersect(const Ray3f &ray_, Intersection &its, bool shadowRay) const
{
    if(!boundingBox.rayIntersect(ray_)) return false;
    if(isLeaf)
    {
        bool foundIntersection = false;  // Was an intersection found so far?
        uint32_t f = (uint32_t) -1;      // Triangle index of the closest intersection
        
        Ray3f ray(ray_); /// Make a copy of the ray (we will need to update its '.maxt' value)
        
        /* Brute force search through all triangles */
        for (uint32_t idx = 0; idx < triangleIndices.size(); ++idx) {
            float u, v, t;
            if (pMesh->rayIntersect(triangleIndices[idx], ray, u, v, t)) {
                /* An intersection was found! Can terminate
                 immediately if this is a shadow ray query */
                if (shadowRay)
                    return true;
                if(t < its.t)
                {
                    ray.maxt = its.t = t;
                    its.uv = Point2f(u, v);
                    its.mesh = pMesh;
                    f = triangleIndices[idx];
                    foundIntersection = true;
                }
            }
        }
        
        if (foundIntersection) {
            /* At this point, we now know that there is an intersection,
             and we know the triangle index of the closest such intersection.
             
             The following computes a number of additional properties which
             characterize the intersection (normals, texture coordinates, etc..)
             */
            
            /* Find the barycentric coordinates */
            Vector3f bary;
            bary << 1-its.uv.sum(), its.uv;
            
            /* References to all relevant mesh buffers */
            const Mesh *mesh   = its.mesh;
            const MatrixXf &V  = mesh->getVertexPositions();
            const MatrixXf &N  = mesh->getVertexNormals();
            const MatrixXf &UV = mesh->getVertexTexCoords();
            const MatrixXu &F  = mesh->getIndices();
            
            /* Vertex indices of the triangle */
            uint32_t idx0 = F(0, f), idx1 = F(1, f), idx2 = F(2, f);
            
            Point3f p0 = V.col(idx0), p1 = V.col(idx1), p2 = V.col(idx2);
            
            /* Compute the intersection positon accurately
             using barycentric coordinates */
            its.p = bary.x() * p0 + bary.y() * p1 + bary.z() * p2;
            
            /* Compute proper texture coordinates if provided by the mesh */
            if (UV.size() > 0)
                its.uv = bary.x() * UV.col(idx0) +
                bary.y() * UV.col(idx1) +
                bary.z() * UV.col(idx2);
            
            /* Compute the geometry frame */
            its.geoFrame = Frame((p1-p0).cross(p2-p0).normalized());
            
            if (N.size() > 0) {
                /* Compute the shading frame. Note that for simplicity,
                 the current implementation doesn't attempt to provide
                 tangents that are continuous across the surface. That
                 means that this code will need to be modified to be able
                 use anisotropic BRDFs, which need tangent continuity */
                
                its.shFrame = Frame(
                                    (bary.x() * N.col(idx0) +
                                     bary.y() * N.col(idx1) +
                                     bary.z() * N.col(idx2)).normalized());
            } else {
                its.shFrame = its.geoFrame;
            }
        }
        
        return foundIntersection;
    }
    else
    {
        bool foundIntersection = false;
        OcTreeNode* sorted_children[8];
        for(int index = 0; index < 8; ++index) sorted_children[index] = children[index];
        for(int index = 0; index < 8; ++index)
        {
            std::stable_sort(sorted_children, sorted_children+8, [ray_](const OcTreeNode* first, const OcTreeNode* second) {
                float firstNearT(MAXFLOAT), secondNearT(MAXFLOAT);
                float firstFarT(MAXFLOAT), secondFarT(MAXFLOAT);
                first->GetBoundingBox().rayIntersect(ray_, firstNearT, firstFarT);
                second->GetBoundingBox().rayIntersect(ray_, secondNearT, secondFarT);
                return firstNearT < secondNearT;
            });
        }
        for(int index = 0; index < 8; ++index)
        {
            foundIntersection = sorted_children[index]->rayIntersect(ray_, its, shadowRay) || foundIntersection;
            if(foundIntersection) break;
        }
        return foundIntersection;
    }
}

void OcTreeNode::SetMeshTriangles(const Mesh* mesh, const std::vector<int>& triangles)
{
    pMesh = mesh;
    triangleIndices = triangles;
}

NORI_NAMESPACE_END
