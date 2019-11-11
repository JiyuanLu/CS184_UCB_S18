#include "student_code.h"
#include "mutablePriorityQueue.h"


using namespace std;

namespace CGL
{
  Vector2D lerp(Vector2D x, Vector2D y, float t){
	return (1 - t) * x + t * y;
  }

  Vector3D lerp3D(Vector3D x, Vector3D y, double t){
	return (1 - t) * x + t * y;
  }

  void BezierCurve::evaluateStep(){
  {
    std::vector<Vector2D> newLevel;
    if(evaluatedLevels.size() == 1){
	for(int i = 0; i < controlPoints.size() -1 ; ++i){
		newLevel.push_back(lerp(controlPoints[i],controlPoints[i+1],t));
	}
    }
    else{
	for(int i = 0; i < evaluatedLevels.back().size() - 1; ++i){
		newLevel.push_back(lerp(evaluatedLevels.back()[i], evaluatedLevels.back()[i+1], t));
	}
    }
    evaluatedLevels.push_back(newLevel);
}
    // TODO Part 1.
    // Perform one step of the Bezier curve's evaluation at t using de Casteljau's algorithm for subdivision.
    // Store all of the intermediate control points into the 2D vector evaluatedLevels.
  }


  Vector3D BezierPatch::evaluate(double u, double v) const
  {
    // TODO Part 2.
    // Evaluate the Bezier surface at parameters (u, v) through 2D de Casteljau subdivision.
    // (i.e. Unlike Part 1 where we performed one subdivision level per call to evaluateStep, this function
    // should apply de Casteljau's algorithm until it computes the final, evaluated point on the surface)
    std::vector<Vector3D> p;
    for(int i = 0; i < controlPoints.size(); ++i){
	p.push_back(evaluate1D(controlPoints[i], u));
    }
    return evaluate1D(p,v);
  }

  Vector3D BezierPatch::evaluate1D(std::vector<Vector3D> points, double t) const
  {
    // TODO Part 2.
    // Optional helper function that you might find useful to implement as an abstraction when implementing BezierPatch::evaluate.
    // Given an array of 4 points that lie on a single curve, evaluates the Bezier curve at parameter t using 1D de Casteljau subdivision.
    std::vector<std::vector<Vector3D>> evaluatedLevels;
    evaluatedLevels.push_back(points);
    while(evaluatedLevels.size() != points.size()){
	std::vector<Vector3D> newLevel;
    	for(int i = 0; i < evaluatedLevels.back().size(); ++i){
		newLevel.push_back(lerp3D(evaluatedLevels.back()[i], evaluatedLevels.back()[i+1],t));
   	 }
    	evaluatedLevels.push_back(newLevel);
    }
    return evaluatedLevels.back()[0];
  }



  Vector3D Vertex::normal( void ) const
  {
    // TODO Part 3.
    // TODO Returns an approximate unit normal at this vertex, computed by
    // TODO taking the area-weighted average of the normals of neighboring
    // TODO triangles, then normalizing.
    Vector3D N(0., 0., 0.);
    HalfedgeCIter h = halfedge();
    h = h -> twin();
    HalfedgeCIter h_orig = h;
    do{
	Vector3D e1 = h->next()->vertex()->position - h->vertex()->position;
	h = h -> next();
	Vector3D e2 = h->next()->vertex()->position - h->vertex()->position;
	N += cross(e1, e2);
        h = h -> twin();
    }while(h != h_orig);	
    return N.unit();
  }

  EdgeIter HalfedgeMesh::flipEdge( EdgeIter e0 )
  {
    // TODO Part 4.
    // TODO This method should flip the given edge and return an iterator to the flipped edge.
    if(!e0->halfedge()->face()->isBoundary() && !e0->halfedge()->twin()->face()->isBoundary()){
	HalfedgeIter h0 = e0->halfedge();
	HalfedgeIter h1 = h0->next();
	HalfedgeIter h2 = h1->next();
	HalfedgeIter h3 = h0->twin();
	HalfedgeIter h4 = h3->next();
	HalfedgeIter h5 = h4->next();
	HalfedgeIter h6 = h1->twin();
	HalfedgeIter h7 = h2->twin();
	HalfedgeIter h8 = h4->twin();
	HalfedgeIter h9 = h5->twin();

	EdgeIter e0 = h0->edge();
	EdgeIter e1 = h1->edge();
	EdgeIter e2 = h2->edge();
	EdgeIter e3 = h4->edge();
	EdgeIter e4 = h5->edge();

	VertexIter v0 = h0->vertex();
	VertexIter v1 = h3->vertex();
	VertexIter v2 = h2->vertex();
	VertexIter v3 = h5->vertex();

	FaceIter f0 = h0->face();
	FaceIter f1 = h3->face();

	// flip
	h0->next() = h1;
	h0->twin() = h3;
	h0->edge() = e0;
	h0->vertex() = v3;
	h0->face() = f0;

	h1->next() = h2;
	h1->twin() = h7;
	h1->edge() = e2;
	h1->vertex() = v2;
	h1->face() = f0;

	h2->next() = h0;
	h2->twin() = h8;
	h2->edge() = e3;
	h2->vertex() = v0;
	h2->face() = f0;

	h3->next() = h4;
	h3->twin() = h0;
	h3->edge() = e0;
	h3->vertex() = v2;
	h3->face() = f1;

	h4->next() = h5;
	h4->twin() = h9;
	h4->edge() = e4;
	h4->vertex() = v3;
	h4->face() = f1;

	h5->next() = h3;
	h5->twin() = h6;
	h5->edge() = e1;
	h5->vertex() = v1;
	h5->face() = f1;
	
	h6->twin() = h5;
	h7->twin() = h1;
	h8->twin() = h2;
	h9->twin() = h4;

	e0->halfedge() = h0;
	e1->halfedge() = h5;
	e2->halfedge() = h1;
	e3->halfedge() = h2;
	e4->halfedge() = h4;
	
	v0->halfedge() = h2;
	v1->halfedge() = h5;
	v2->halfedge() = h3;
	v3->halfedge() = h0;

	f0->halfedge() = h0;
	f1->halfedge() = h3;
    }
    return e0;
  }

  VertexIter HalfedgeMesh::splitEdge( EdgeIter e0 )
  {
    // TODO Part 5.
    // TODO This method should split the given edge and return an iterator to the newly inserted vertex.
    // TODO The halfedge of this vertex should point along the edge that was split, rather than the new edges.
    if(!e0->halfedge()->face()->isBoundary() && !e0->halfedge()->twin()->face()->isBoundary()){
	// collect pointers, similar to part 4
	// the 4 outside halfedges are not needed
	HalfedgeIter h0 = e0->halfedge();
	HalfedgeIter h1 = h0->next();
	HalfedgeIter h2 = h1->next();
	HalfedgeIter h3 = h0->twin();
	HalfedgeIter h4 = h3->next();
	HalfedgeIter h5 = h4->next();
	HalfedgeIter h6 = h1->twin();
      	HalfedgeIter h7 = h2->twin();
      	HalfedgeIter h8 = h4->twin();
      	HalfedgeIter h9 = h5->twin();

	EdgeIter e0 = h0->edge();
	EdgeIter e1 = h1->edge();
	EdgeIter e2 = h2->edge();
	EdgeIter e3 = h4->edge();
	EdgeIter e4 = h5->edge();

	VertexIter v0 = h0->vertex();
	VertexIter v1 = h3->vertex();
	VertexIter v2 = h2->vertex();
	VertexIter v3 = h5->vertex();

	FaceIter f0 = h0->face();
	FaceIter f1 = h3->face();
	
	// allocate new pointers
	Vector3D x = e0->halfedge()->vertex()->position;
	Vector3D y = e0->halfedge()->twin()->vertex()->position;
	Vector3D mid = (x+y)/2;
	VertexIter m = newVertex();
	m->position = mid;
	m->isNew = true;
	
	EdgeIter e5 = newEdge();
	EdgeIter e6 = newEdge();
	EdgeIter e7 = newEdge();

	e6->isNew = true;
	e7->isNew = true;

	HalfedgeIter h10 = newHalfedge();
	HalfedgeIter h11 = newHalfedge();
	HalfedgeIter h12 = newHalfedge();
	HalfedgeIter h13 = newHalfedge();
	HalfedgeIter h14 = newHalfedge();
	HalfedgeIter h15 = newHalfedge();

	FaceIter f2 = newFace();
	FaceIter f3 = newFace();

	// reassign
 	h0->next() = h13;
     	h0->twin() = h10;
      	h0->vertex() = v0;
      	h0->edge() = e0;
      	h0->face() = f0;

      	h1->next() = h12;
	h1->twin() = h6;
      	h1->vertex() = v1;
      	h1->edge() = e1;
      	h1->face() = f2;

      	h2->next() = h0;
	h2->twin() = h7;
      	h2->vertex() = v2;
      	h2->edge() = e2;
      	h2->face() = f0;

      	h3->next() = h15;
      	h3->twin() = h11;
      	h3->vertex() = v1;
      	h3->edge() = e5;
      	h3->face() = f1;

      	h4->next() = h14;
	h4->twin() = h8;
      	h4->vertex() = v0;
      	h4->edge() = e3;
      	h4->face() = f3;

      	h5->next() = h3;
	h5->twin() = h9;
      	h5->vertex() = v3;
      	h5->edge() = e4;
      	h5->face() = f1;

	h6->next() = h6->next();
        h6->twin() = h1;
        h6->vertex() = v2;
        h6->edge() = e1;
        h6->face() = h6->face();

        h7->next() = h7->next();
        h7->twin() = h2;
        h7->vertex() = v0;
        h7->edge() = e2;
        h7->face() = h7->face();

        h8->next() = h8->next();
        h8->twin() = h4;
        h8->vertex() = v3;
        h8->edge() = e3;
        h8->face() = h8->face();

        h9->next() = h9->next();
        h9->twin() = h5;
        h9->vertex() = v1;
        h9->edge() = e4;
        h9->face() = h9->face();

      	h10->next() = h4;
      	h10->twin() = h0;
      	h10->vertex() = m;
      	h10->edge() = e0;
      	h10->face() = f3;

      	h11->next() = h1;
	h11->twin() = h3;
      	h11->vertex() = m;
      	h11->edge() = e5;
      	h11->face() = f2;

      	h12->next() = h11;
      	h12->twin() = h13;
      	h12->vertex() = v2;
      	h12->edge() = e7;
      	h12->face() = f2;

      	h13->next() = h2;
      	h13->twin() = h12;
      	h13->vertex() = m;
      	h13->edge() = e7;
      	h13->face() = f0;

      	h14->next() = h10;
      	h14->twin() = h15;
      	h14->vertex() = v3;
      	h14->edge() = e6;
      	h14->face() = f3;

      	h15->next() = h5;
      	h15->twin() = h14;
      	h15->vertex() = m;
      	h15->edge() = e6;
      	h15->face() = f1;

      	m->halfedge() = h10;  
      	v0->halfedge() = h0;
      	v1->halfedge() = h3;
      	v2->halfedge() = h12;
      	v3->halfedge() = h14;

      
      	e0->halfedge() = h0;
      	e1->halfedge() = h1;
      	e2->halfedge() = h2;
      	e3->halfedge() = h4;
      	e4->halfedge() = h5;
      	e5->halfedge() = h3;
      	e6->halfedge() = h14;
      	e7->halfedge() = h12;

      	f0->halfedge() = h0;
      	f1->halfedge() = h3;
      	f2->halfedge() = h11;
      	f3->halfedge() = h10;

      	return m;
    }
  }



  void MeshResampler::upsample( HalfedgeMesh& mesh )
  {
    // TODO Part 6.
    // This routine should increase the number of triangles in the mesh using Loop subdivision.
    // Each vertex and edge of the original surface can be associated with a vertex in the new (subdivided) surface.
    // Therefore, our strategy for computing the subdivided vertex locations is to *first* compute the new positions
    // using the connectity of the original (coarse) mesh; navigating this mesh will be much easier than navigating
    // the new subdivided (fine) mesh, which has more elements to traverse. We will then assign vertex positions in
    // the new mesh based on the values we computed for the original mesh.


    // TODO Compute new positions for all the vertices in the input mesh, using the Loop subdivision rule,
    for(VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); ++v){
	Vector3D original_position = v->position;
	Vector3D neighboring_positions(0,0,0);
	HalfedgeCIter h = v->halfedge();
	int num = 0;	//count # of neighbors
	do{
		neighboring_positions += h->twin()->vertex()->position;
		h = h->twin()->next();
		++num;
	}while(h != v->halfedge());
	double u = 0.0;
	if(num == 3)
		u = 3.0/16.0;
	else
		u = 3.0/(8.0*num);
	v->newPosition = (1-num*u) * original_position + u * neighboring_positions;
	v->isNew = false;
    }
		
    // TODO and store them in Vertex::newPosition. At this point, we also want to mark each vertex as being
    // TODO a vertex of the original mesh.


    // TODO Next, compute the updated vertex positions associated with edges, and store it in Edge::newPosition.
	for(EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); ++e){
		Vector3D AB_positions(0,0,0);
		Vector3D CD_positions(0,0,0);
		AB_positions = e->halfedge()->vertex()->position + e->halfedge()->twin()->vertex()->position;
		CD_positions = e->halfedge()->next()->twin()->vertex()->position + e->halfedge()->twin()->next()->twin()->vertex()->position;
		e->newPosition = 3.0/8.0 * AB_positions + 1.0/8.0 * CD_positions;
		e->isNew = false;
	}
    // TODO Next, we're going to split every edge in the mesh, in any order.  For future
    // TODO reference, we're also going to store some information about which subdivided
    // TODO edges come from splitting an edge in the original mesh, and which edges are new,
    // TODO by setting the flat Edge::isNew.  Note that in this loop, we only want to iterate
    // TODO over edges of the original mesh---otherwise, we'll end up splitting edges that we
    // TODO just split (and the loop will never end!)
    	
	for(EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd();){
		EdgeIter nextEdge = e;
		++nextEdge;
		if(!e->halfedge()->vertex()->isNew && !e->halfedge()->twin()->vertex()->isNew){
			VertexIter newVertex = mesh.splitEdge(e);
			newVertex->newPosition = e->newPosition;
		}
		e = nextEdge;
	}

    // TODO Now flip any new edge that connects an old and new vertex.
	for(EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd();){
		EdgeIter nextEdge = e;
		++nextEdge;
		// did I set e->isNew = true for new edge e when splitting?
		if(e->isNew){
			bool b1 = e->halfedge()->vertex()->isNew;
			bool b2 = e->halfedge()->twin()->vertex()->isNew;
			if((b1 && !b2) || (!b1 && b2)){
				mesh.flipEdge(e);
			}
		}
		e = nextEdge;
	}
			
    // TODO Finally, copy the new vertex positions into final Vertex::position.
	for(VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); ++v){
		v->position = v->newPosition;
	}
  }
	
}
