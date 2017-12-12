//
//#ifndef ALGORITHM_H
//#define ALGORITHM_H
//
//#include <common/meshmodel.h>
//#include <common/interfaces.h>
//#include "sparse_matrix.h"
//#include "pcg_solver.h"
//#include <vector>
//#include <math.h>
//using namespace std;
//
//namespace spin{
//
//void  buildB(MeshModel &m,std::map<int,vcg::Point3f>& map,std::vector<double>& b){
//	 // visit each face
//   for(int i = 0; i < m.cm.face.size(); i++) {
//	   CFaceO &face = m.cm.face[i];	   
//	   int col = face.Index();
//	   vcg::Point3f& kexi = map[face.Index()];
//        for(int j = 0;j < 3; j++) {
//			// get vertex positions
//			CVertexO* v0 = face.V( (j) % 3 );
//			CVertexO* v1 = face.V( (j+1) % 3 );
//			CVertexO* v2 = face.V( (j+2) % 3 );
//
//			// get vertex indices
//			int row = v0->Index();
//			// int k1 = v1->Index();
//			// int k2 = v2->Index();
//
//			// compute cotangent of the angle at the current vertex
//			// (equal to cosine over sine, which equals the dot
//			// product over the norm of the cross product)
//			vcg::Point3f u12 = v2->P() - v1->P();
//			vcg::Point3f u01 = v1->P() - v0->P();
//			vcg::Point3f normal_face=(u01^u12).Normalize();
//				 
//			vcg::Point3f nabla_phi0 = u12^normal_face;
//		// spin::Quaternion q(0,nabla_phi0.X()/2,nabla_phi0.Y()/2,nabla_phi0.Z()/2);
//			// add contribution of this cotangent to the matrix
//			b[row] += nabla_phi0*kexi;
//		}
//	}
//	
//}
//
//
//void  buildLaplacian(MeshModel &m,SparseMatrixd&L)
//// builds the cotan-Laplace operator
//{
//   // allocate a sparse |V|x|V| matrix
//  
//	int nV = m.cm.VN();
//	int nF = m.cm.FN();
//    L.resize(nV);
//
//   // visit each face
//   for(int i = 0; i < m.cm.face.size(); i++) {
//	   CFaceO &face = m.cm.face[i];	   
//	   int col = face.Index();
//        for(int j = 0;j < 3; j++) {
//			// get vertex positions
//			CVertexO* v0 = face.V( (j) % 3 );
//			CVertexO* v1 = face.V( (j+1) % 3 );
//			CVertexO* v2 = face.V( (j+2) % 3 );
//
//			// get vertex indices
//			int row = v0->Index();
//			int k1 = v1->Index();
//			int k2 = v2->Index();
//
//			// compute cotangent of the angle at the current vertex
//			// (equal to cosine over sine, which equals the dot
//			// product over the norm of the cross product)       
//			vcg::Point3f u1 = v1->P() - v0->P();
//			vcg::Point3f u2 = v2->P() - v0->P();
//			double cotAlpha = (u1*u2)/(u1^u2).Norm();
//
//			// add contribution of this cotangent to the matrix
//		/* L( k1, k2 ) -= cotAlpha / 2.;
//			L( k2, k1 ) -= cotAlpha / 2.;
//			L( k1, k1 ) += cotAlpha / 2.;
//			L( k2, k2 ) += cotAlpha / 2.;*/
//
//			L.add_to_element(k1,k2,-cotAlpha / 2.);
//			L.add_to_element(k2,k1,-cotAlpha / 2.);
//			L.add_to_element(k1,k2,cotAlpha / 2.);
//			L.add_to_element(k2,k1,-cotAlpha / 2.);
//		 }
//	}
//}
//
//
//void solve( SparseMatrixd& A,
//                     vector<double>& x,
//                     vector<double>& b,
//                     bool precondition )
//// solves the linear system Ax = b where A is positive-semidefinite
//{
//   if( precondition == false )
//   {
//      cerr << "WARNING: using basic CG solver with diagonal preconditioner -- may be (very) slow!" << endl;
//   }
//
//   int n = x.size();
//   vector<double> result( n );
//   //vector<float> rhs( n );
//
//   // setup solver
//   PCGSolver<double> pcg;
//   const int max_iterations = 2000;
//   //const int max_iterations = 200;
//   const double tolerance_factor = 1e-7;
//   const double mic_parameter = .97;
//   const double min_diagonal_ratio = .25;
//   double residual;
//   int iterations;
//
//   pcg.set_solver_parameters( tolerance_factor, max_iterations, mic_parameter, min_diagonal_ratio, precondition );
//   
//   // solve real linear system
//   pcg.solve( A, b, result, residual, iterations );
//
//   cout << "Linear solver achieved a residual of " << residual;
//   cout << " after " << iterations << " iterations." << endl;
//}
//
//void countCurlfree(MeshModel &m,std::map<int,vcg::Point3f>& kexi,vector<double>& x){
//	SparseMatrixd L;// 系数矩阵|V|*|F|，为（4）中左侧的部分
////    vector<float> x;// 解向量共|F|个分量
//	x.resize(m.cm.VN());
//	std::vector<double> b;
//	b.resize(m.cm.VN());
//	buildB(m,kexi,b);
//	buildLaplacian(m,L);
//	solve( L, x,b ,false);
//
//}
//
//}// end of namespace
//#endif