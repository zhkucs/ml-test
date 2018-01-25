#ifndef _ALG_BASE_H
#define _ALG_BASE_H

#include <eigen/Sparse>
#include <common/meshmodel.h>
#include <common/interfaces.h>
#include <vector>
#include <math.h>
#include "wrap/gl/gl_geometry.h"

using namespace std;
using namespace vcg;
//typedef Eigen::Triplet<double> T;
//typedef Eigen::SparseMatrix<double> SparseMatrixType;
#define EPSLON 0.001
#define PI 3.141592653589
#define rad2angle( r) r*180/PI;
# define KH_DIFF_THRESHOLD 0.002
#define  equal_by_threshold(v1,v2) (abs(v1-v2) < KH_DIFF_THRESHOLD)

struct Color{
		Color(float r,float g,float b){
			_r =r;
			_g = g;
			_b = b;
		}
		float _r,_g,_b;

	};
struct ShowParameter{
		ShowParameter(){}
		ShowParameter(double r,double d,double s){
			_r = r;
			_dotSize = d;
			_singularityDot = s;
		}
		double _r;// 三角形的外心半径，可以用来控制向量长度
		double _dotSize;// 针头半径大小
		double _singularityDot;// 奇异点半径
	};
namespace gdut_base{
	//函数对象 用于对pair进行hash

	enum VertexType{
		SOURCE,
		SINK,
		TRIVAL
	};
	struct pairhash{
		template<class T1, class T2>
		size_t operator()(const pair<T1, T2> &x) const{ //为什么去掉最后的const编译无法通过？？？？？？？
			hash<T1> h1;
			hash<T2> h2;
			return h1(x.first) ^ h2(x.second); 
		}
	};
	

	

	inline float countArea(CFaceO& f){
		vcg::Point3f p_i = f.P(0);
		vcg::Point3f p_j = f.P(1);
		vcg::Point3f p_k = f.P(2);
		vcg::Point3f vij=p_j-p_i;  
		vcg::Point3f vik=-p_k-p_i;
		return (vij ^ vik).Norm()/2;
	}

	inline void countNablaOfFace(CFaceO&f,double s0,double s1,double s2,vcg::Point3f&result){
		vcg::Point3f p_i = f.P0(0); 
		vcg::Point3f p_j = f.P0(1); 
		vcg::Point3f p_k = f.P0(2); 

		//梯度的计算参考了 http://blog.csdn.net/zdy0_2004/article/details/49615919
		double double_area =2*countArea(f);

		vcg::Point3f normalf = NormalizedNormal<CFaceO>(f); 

		vcg::Point3f phi_j = normalf^(p_i-p_k)/double_area;
		vcg::Point3f phi_k = normalf^(p_j-p_i)/double_area;

		result=phi_j*(s1-s0)+phi_k*(s2-s0);	
	}

	inline void buildKexi(MeshModel &m,std::map<int,Point3f>& kexi){
		for(CMeshO::FaceIterator fi=m.cm.face.begin();fi!=m.cm.face.end(); ++fi)	{
			CFaceO f = *fi;  	
			f.V0(0)->Base().CurvatureEnabled=true;


			//	float areaf=0.5f*normalf.Norm();  
			//	vcg::Point3f nr = normalf.Normalize();
			vcg::Point3f bc = Barycenter(f);
			float kh_i = f.V0(0)->Kh();
			float kh_j = f.V0(1)->Kh();
			float kh_k = f.V0(2)->Kh();
			//double r = Distance(f.P0(0),bc);

			vcg::Point3f nabla_f;
			countNablaOfFace(f,kh_i,kh_j,kh_k,nabla_f);

			vcg::Point3f start = Barycenter(f);
			vcg::Point3f end =  start + (nabla_f*4);
			//vcg::Point3f newEnd = standardize(start,end,r);// 画图为了好看，将向量缩放到三角形范围内，实际梯度的计算仍用回start-->end

			kexi.insert(pair<int,Point3f&>((*fi).Index(),end-start));
			//start
		}
	}

	// 根据三点计算在pi处的拉格朗日基函数的梯度；pi,pj,pk的次序为逆时针方向。
	inline void countPhi(vcg::Point3f& p_i,vcg::Point3f& p_j,vcg::Point3f& p_k,vcg::Point3f&result){
		vcg::Point3f vij=p_j-p_i;  
		vcg::Point3f vik=-p_k-p_i;
		vcg::Point3f normal = (vik^vij);
		result= (p_j-p_k)^normal.Normalize()/normal.Norm();
	}

	//遍历的算法见 http://vcg.isti.cnr.it/vcglib/adjacency.html
	inline void extremum_kh_1_ring(CVertexO * v,std::set<int>& v_Chected,pair<int,int>& max_min_indice,pair<float,float>& max_min_kh){
//		float kh_v = v->Kh();		
		CMeshO::FacePointer fp = v->VFp();
		CFaceO* start = &fp[0];
		vcg::face::Pos<CFaceO> pos(start,v);// constructor that takes face, edge and vertex
		pos.FlipV();// 得到共边&共面的另一个v
		do
		{
			if(v_Chected.count(pos.v->Index())){
				break;
			}

			float current_kh = pos.v->Kh();	
			if(current_kh > max_min_kh.first) 
				max_min_indice.first = pos.v->Index();
			else if(current_kh < max_min_kh.second) 
				max_min_indice.second = pos.v->Index();
		
			//pos.FlipV();// 变回来
			pos.FlipF();// 得到共边，共点的下一个cell（面不同）
			pos.FlipE();// 得到共面，共点的下一个cell（边不同）

		}while(pos.f!=start);
	}

	inline VertexType isSupport(CVertexO * v,int kernelIndex,float kh){
		CMeshO::FacePointer fp = v->VFp();
		CFaceO* start = &fp[0];
		vcg::face::Pos<CFaceO> pos(start,v);// constructor that takes face, edge and vertex
		pos.FlipV();// 得到共边&共面的另一个v
		float kh_0 = pos.v->Kh();
		pos.FlipV();// 变回来
		if(kh > kh_0){// 核心大于周边
			do
			{	
				// 下一个cell
				pos.FlipF();// 得到共边，共点的下一个cell（面不同）
				pos.FlipE();// 得到共面，共点的下一个cell（边不同）
				pos.FlipV();

				if(kernelIndex != pos.v->Index()){// pos.v不是核心才处理
					float current_kh = pos.v->Kh();	
					if(kh<current_kh || equal_by_threshold(kh,current_kh))
						return TRIVAL;
				}

				pos.FlipV();
			}while(pos.f!=start);
			return SOURCE;
		}else{// 核心小于周边
			do
			{		
				// 下一个cell
				pos.FlipF();// 得到共边，共点的下一个cell（面不同）
				pos.FlipE();// 得到共面，共点的下一个cell（边不同）
				pos.FlipV();

				if(kernelIndex != pos.v->Index()){// pos.v不是核心才处理
					float current_kh = pos.v->Kh();	
					if(kh>current_kh|| equal_by_threshold(kh,current_kh))
						return TRIVAL;
				}
		
				pos.FlipV();

			}while(pos.f!=start);
			return SINK;
		}
	}

	inline VertexType extremum_kh_2_ring(CVertexO * v){
		float kh_v = v->Kh();		
		VertexType vt0 = isSupport(v,v->Index(),kh_v);
		if(vt0 == TRIVAL) return TRIVAL;

		// 取得第一个面开始遍历
		CMeshO::FacePointer fp = v->VFp();
		CFaceO* start = &fp[0];
		vcg::face::Pos<CFaceO> pos(start,v);// constructor that takes face, edge and vertex
		//pos.FlipV();// 得到共边&共面的另一个v
		do
		{		
			// 下一个cell
			pos.FlipF();// 得到共边，共点的下一个cell（面不同）
			pos.FlipE();// 得到共面，共点的下一个cell（边不同）
			pos.FlipV();

			VertexType vt1 = isSupport(pos.v,v->Index(),kh_v);
			if(vt1 == TRIVAL) 
				return TRIVAL;
		
			pos.FlipV();

		}while(pos.f!=start);

		return vt0;
	}
	inline void drawArrow(vcg::Point3f &origin,vcg::Point3f &dst,Color &color)
	{//////////////////////////////////////////
		vcg::Point3f z(0,0,1);
		vcg::Point3f x_axis(1,0,0);
		vcg::Point3f y(0,1,0);	
		vcg::Point3f zero(0,0,0);

		vcg::Point3f diff(dst-origin);

		vcg::Point3f diff_xoy(diff.X(),diff.Y(),0);

		vcg::Point3f axis2 =diff_xoy^diff;

		// 在xoy面上画箭头两翼（尖端为dst）
		double rad = 0.15;	
		vcg::Point3f g = dst - origin;
		float wing = g.Norm()*0.3;



		// 两次旋转角度
//		double angle1 = rad2angle(acos(diff.X()/diff_xoy.Norm()));// diff_xoy 与x轴夹角	
//		double angle2 = rad2angle(asin(diff.Z()/diff.Norm()));// diff与diff_xoy的夹角

		vcg::Point3f rotateAxis = x_axis ^ diff;
		double r = asin(rotateAxis.Norm()/diff.Norm());
		double angle3 = rad2angle(r);// diff与diff_xoy的夹角

		/// 画带箭头的直线
		glLineWidth(1.0);	
		glEnable(GL_LINE_SMOOTH);
		/*glBegin(GL_LINES);
		glVertex(origin);
		glVertex(dst);	
		glEnd();*/
		float len = g.Norm();
		glPushMatrix();			
		glTranslatef(origin.X(), origin.Y(), origin.Z());
		//glRotatef(angle1,z.X(),z.Y(),z.Z());
		//glRotatef(angle2,axis2.X(),axis2.Y(),axis2.Z());

		glRotatef(angle3,rotateAxis.X(),rotateAxis.Y(),rotateAxis.Z());

		// 画箭头干 
		glBegin(GL_LINES);
		glColor3f(color._r,color._g,color._b);
		glVertex(zero);
		glVertex(vcg::Point3f(len,0,0));	
		glEnd();

		vcg::Point3f rear1(len-wing*cos(rad),wing*sin(rad),0);
		vcg::Point3f rear2(len-wing*cos(rad),-wing*sin(rad),0);
		glBegin(GL_LINE_STRIP);		
		//glColor3f(0.0,1.0,0);
		glVertex(rear1);
		glVertex(vcg::Point3f(len,0,0));
		//glColor3f(0.0,0,1.0);
		glVertex(rear2);
		glEnd();
		glPopMatrix();
		glFlush();                                                                                                             
		glDisable(GL_LINE_SMOOTH);
	}

	inline void drawStick(vcg::Point3f &origin,vcg::Point3f &dst,Color &color,ShowParameter& para)
	{//////////////////////////////////////////
		// 画箭头干
		glBegin(GL_LINES);
		glColor3f(color._r,color._g,color._b);
		glVertex(origin);
		glVertex(dst);	
		glEnd();

		GLdouble r = para._dotSize;// bunny模型要用0.05，kitty要用0.0006
		GLint m = 8;
		GLint n = 8;
		glPushMatrix();			
		glTranslatef(origin.X(), origin.Y(), origin.Z());
		::glutSolidSphere(r, m, n); 
		glPopMatrix();

		//glEnd();
		glFlush();                                                                                                             
		glDisable(GL_LINE_SMOOTH);
	}

	inline void drawStickMapOnface(vcg::Point3f &origin,vcg::Point3f &dst,vcg::Point3f &normal,Color &color)
	{//////////////////////////////////////////
		// 画箭头干
		glBegin(GL_LINES);
		glColor3f(color._r,color._g,color._b);
		glVertex(origin);
		glVertex(dst);	
		glEnd();

		GLdouble r = 0.05;
		GLint m = 5;
		GLint n = 5;
		glPushMatrix();			
		glTranslatef(origin.X(), origin.Y(), origin.Z());
		::glutSolidSphere(r, m, n); 
		glPopMatrix();

		glEnd();
		glFlush();                                                                                                             
		glDisable(GL_LINE_SMOOTH);
	}

	inline void drawPoints(std::set<vcg::Point3f> &pts,Color &color,ShowParameter& para)
	{//////////////////////////////////////////

		GLdouble r = para._singularityDot;
		GLint m = 5;
		GLint n = 5;
		glColor3f(color._r,color._g,color._b);
		for(set<vcg::Point3f>::iterator it = pts.begin();it!= pts.end();it++){
			glPushMatrix();			
			glTranslatef((*it).X(), (*it).Y(), (*it).Z());
			::glutSolidSphere(r, m, n); 
			glPopMatrix();	
		}
		glFlush();                                                                                                             
		glDisable(GL_LINE_SMOOTH);
	}

	inline void drawArrowOnFace(vcg::Point3f &origin,vcg::Point3f &dst,vcg::Point3f &normal,Color &color)
	{//////////////////////////////////////////
		vcg::Point3f z(0,0,1);
		vcg::Point3f x_axis(1,0,0);
		vcg::Point3f y(0,1,0);	
		vcg::Point3f zero(0,0,0);

		// 在xoy面上画箭头两翼（尖端为dst）
		double rad = 0.15;	
		vcg::Point3f g = dst - origin;
		float wing = g.Norm()*0.3;


		// 两个旋转轴
		vcg::Point3f axis_n_z = z ^ normal;
		vcg::Point3f axis_g_x = z;
		// 两次旋转角度
		double angle_n_z = rad2angle(acos(normal*z/normal.Norm()));// n与z的夹角	
		double angle_g_x = rad2angle(acos((g^normal)*y/g.Norm()));// g与x的夹角

		/// 画带箭头的直线
		glLineWidth(1.0);	
		glEnable(GL_LINE_SMOOTH);
		/*glBegin(GL_LINES);
		glVertex(origin);
		glVertex(dst);	
		glEnd();*/
		float len = g.Norm();
		glPushMatrix();			
		glTranslatef(origin.X(), origin.Y(), origin.Z());
		glRotatef(angle_n_z,axis_n_z.X(),axis_n_z.Y(),axis_n_z.Z());
		glRotatef(angle_g_x,axis_g_x.X(),axis_g_x.Y(),axis_g_x.Z());

		// 画箭头干
		glBegin(GL_LINES);
		glColor3f(color._r,color._g,color._b);
		glVertex(zero);
		glVertex(vcg::Point3f(len,0,0));	
		glEnd();

		vcg::Point3f rear1(len-wing*cos(rad),wing*sin(rad),0);
		vcg::Point3f rear2(len-wing*cos(rad),-wing*sin(rad),0);
		glBegin(GL_LINE_STRIP);		
		//glColor3f(0.0,1.0,0);
		glVertex(rear1);
		glVertex(vcg::Point3f(len,0,0));
		//glColor3f(0.0,0,1.0);
		glVertex(rear2);
		glEnd();
		glPopMatrix();
		glFlush();                                                                                                             
		glDisable(GL_LINE_SMOOTH);
	}

	// 检查局部曲率最大/小的顶点，并画出
	inline void findSingluarityOnVertex(CVertexO* v,std::set<int>& v_Chected,std::set<Point3f>& s_source,std::set<Point3f>& s_sink){		
			if(v_Chected.count(v->Index())){// 检查顶点是否已经处理过
				return;
			}
			
			//gdut_base::VertexType vt = gdut_base::countCenter(v);
			gdut_base::VertexType vt = gdut_base::extremum_kh_2_ring(v);
			switch(vt){
				case gdut_base::SOURCE:
					//drawStick(v->P() ,v->P(),gdut_base::Red);
					s_source.insert(v->P());
					v_Chected.insert(v->Index());
					return;
				case gdut_base::SINK:
					//drawStick(v->P() ,v->P(),gdut_base::Green);
					s_sink.insert(v->P());
					v_Chected.insert(v->Index());
					return;
			}
	}

	// 根据顶点平均曲率找到2-ring极值顶点
	inline void findSingluarityOnMesh(CMeshO& m,std::set<Point3f>& s_source,std::set<Point3f>& s_sink){
		std::set<int> v_Chected;
		for(CMeshO::VertexIterator vi=m.vert.begin();vi!=m.vert.end(); ++vi){				
			gdut_base::findSingluarityOnVertex(&(*vi),v_Chected,s_source,s_sink);
		}
	}
}

#endif