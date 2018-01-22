#ifndef ALG_BASE_H
#define ALG_BASE_H

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
namespace gdut_base{
	//�������� ���ڶ�pair����hash

	enum VertexType{
		SOURCE,
		SINK,
		TRIVAL
	};
	struct pairhash{
		template<class T1, class T2>
		size_t operator()(const pair<T1, T2> &x) const{ //Ϊʲôȥ������const�����޷�ͨ����������������
			hash<T1> h1;
			hash<T2> h2;
			return h1(x.first) ^ h2(x.second); 
		}
	};
	struct Color{
		Color(float r,float g,float b){
			_r =r;
			_g = g;
			_b = b;
		}
		float _r,_g,_b;

	}Red(1,0,0),Green(0,1,0),Blue(0,0,1),Orange(1,0.5f,0);

	float countArea(CFaceO& f){
		vcg::Point3f p_i = f.P(0);
		vcg::Point3f p_j = f.P(1);
		vcg::Point3f p_k = f.P(2);
		vcg::Point3f vij=p_j-p_i;  
		vcg::Point3f vik=-p_k-p_i;
		return (vij ^ vik).Norm()/2;
	}

	void countNablaOfFace(CFaceO&f,double s0,double s1,double s2,vcg::Point3f&result){
		vcg::Point3f p_i = f.P0(0); 
		vcg::Point3f p_j = f.P0(1); 
		vcg::Point3f p_k = f.P0(2); 

		//�ݶȵļ���ο��� http://blog.csdn.net/zdy0_2004/article/details/49615919
		double double_area =2*countArea(f);

		vcg::Point3f normalf = NormalizedNormal<CFaceO>(f); 

		vcg::Point3f phi_j = normalf^(p_i-p_k)/double_area;
		vcg::Point3f phi_k = normalf^(p_j-p_i)/double_area;

		result=phi_j*(s1-s0)+phi_k*(s2-s0);	
	}

	void buildKexi(MeshModel &m,std::map<int,Point3f>& kexi){
		for(CMeshO::FaceIterator fi=m.cm.face.begin();fi!=m.cm.face.end(); ++fi)	{
			CFaceO f = *fi;  	
			f.V0(0)->Base().CurvatureEnabled=true;


			//	float areaf=0.5f*normalf.Norm();  
			//	vcg::Point3f nr = normalf.Normalize();
			vcg::Point3f bc = Barycenter(f);
			float kh_i = f.V0(0)->Kh();
			float kh_j = f.V0(1)->Kh();
			float kh_k = f.V0(2)->Kh();
			double r = Distance(f.P0(0),bc);

			vcg::Point3f nabla_f;
			countNablaOfFace(f,kh_i,kh_j,kh_k,nabla_f);

			vcg::Point3f start = Barycenter(f);
			vcg::Point3f end =  start + (nabla_f*4);
			//vcg::Point3f newEnd = standardize(start,end,r);// ��ͼΪ�˺ÿ������������ŵ������η�Χ�ڣ�ʵ���ݶȵļ������û�start-->end

			kexi.insert(pair<int,Point3f&>((*fi).Index(),end-start));
			//start
		}
	}

	// �������������pi�����������ջ��������ݶȣ�pi,pj,pk�Ĵ���Ϊ��ʱ�뷽��
	void countPhi(vcg::Point3f& p_i,vcg::Point3f& p_j,vcg::Point3f& p_k,vcg::Point3f&result){
		vcg::Point3f vij=p_j-p_i;  
		vcg::Point3f vik=-p_k-p_i;
		vcg::Point3f normal = (vik^vij);
		result= (p_j-p_k)^normal.Normalize()/normal.Norm();
	}

	//�������㷨�� http://vcg.isti.cnr.it/vcglib/adjacency.html
	void extremum_kh_1_ring(CVertexO * v,pair<int,int>& max_min_indice){
		float kh_v = v->Kh();		
		max_min_indice.first = max_min_indice.second = v->Index();
		CMeshO::FacePointer fp = v->VFp();
		CFaceO* start = &fp[0];
		vcg::face::Pos<CFaceO> pos(start,v);// constructor that takes face, edge and vertex
		do
		{
			pos.FlipV();// �õ�����&�������һ��v

			float current_kh = pos.v->Kh();	
			if(current_kh > kh_v) max_min_indice.first = pos.v->Index();
			else if(current_kh < kh_v) max_min_indice.second = pos.v->Index();
		
			pos.FlipV();// �����
			pos.FlipF();// �õ����ߣ��������һ��cell���治ͬ��
			pos.FlipE();// �õ����棬�������һ��cell���߲�ͬ��

		}while(pos.f!=start);

	}

	void extremum_kh_2_ring(CVertexO * v){
		pair<int,int>& max_min_indice;
		float kh_v = v->Kh();		
		max_min_indice.first = max_min_indice.second = v->Index();
		CMeshO::FacePointer fp = v->VFp();
		CFaceO* start = &fp[0];
		vcg::face::Pos<CFaceO> pos(start,v);// constructor that takes face, edge and vertex
		do
		{
			pos.FlipV();// �õ�����&�������һ��v

			pair<int,int>& mmi(max_min_indice);
			extremum_kh_1_ring(pos.v,mmi);
		
			pos.FlipV();// �����
			pos.FlipF();// �õ����ߣ��������һ��cell���治ͬ��
			pos.FlipE();// �õ����棬�������һ��cell���߲�ͬ��

		}while(pos.f!=start);

	}
	VertexType countCenter(CVertexO * v){
		float kh_v = v->Kh();
		bool all_greater = true;
		bool all_smaller = true;
		CMeshO::FacePointer fp = v->VFp();
		CFaceO* start = &fp[0];
		vcg::face::Pos<CFaceO> pos(start,v);// constructor that takes face, edge and vertex
		do
		{
			pos.FlipV();// �õ�����&�������һ��v

			float current_kh = pos.v->Kh();
			all_greater = (all_greater&&current_kh >= kh_v);
			all_smaller = (all_smaller&&current_kh <= kh_v);
		
			pos.FlipV();// �����
			pos.FlipF();// �õ����ߣ��������һ��cell���治ͬ��
			pos.FlipE();// �õ����棬�������һ��cell���߲�ͬ��

		}while(pos.f!=start);
		if(all_greater) return SINK;
		if(all_smaller) return SOURCE;
		return TRIVAL;

	}






	void drawArrow(vcg::Point3f &origin,vcg::Point3f &dst,Color &color)
	{//////////////////////////////////////////
		vcg::Point3f z(0,0,1);
		vcg::Point3f x_axis(1,0,0);
		vcg::Point3f y(0,1,0);	
		vcg::Point3f zero(0,0,0);

		vcg::Point3f diff(dst-origin);

		vcg::Point3f diff_xoy(diff.X(),diff.Y(),0);

		vcg::Point3f axis2 =diff_xoy^diff;

		// ��xoy���ϻ���ͷ�������Ϊdst��
		double rad = 0.15;	
		vcg::Point3f g = dst - origin;
		float wing = g.Norm()*0.3;



		// ������ת�Ƕ�
		double angle1 = rad2angle(acos(diff.X()/diff_xoy.Norm()));// diff_xoy ��x��н�	
		double angle2 = rad2angle(asin(diff.Z()/diff.Norm()));// diff��diff_xoy�ļн�

		vcg::Point3f rotateAxis = x_axis ^ diff;
		double r = asin(rotateAxis.Norm()/diff.Norm());
		double angle3 = rad2angle(r);// diff��diff_xoy�ļн�

		/// ������ͷ��ֱ��
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

		// ����ͷ�� 
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

	void drawStick(vcg::Point3f &origin,vcg::Point3f &dst,Color &color)
	{//////////////////////////////////////////
		// ����ͷ��
		glBegin(GL_LINES);
		glColor3f(color._r,color._g,color._b);
		glVertex(origin);
		glVertex(dst);	
		glEnd();

		GLdouble r = 0.0006;// bunnyģ��Ҫ��0.05��kittyҪ��0.0006
		GLint m = 8;
		GLint n = 8;
		glPushMatrix();			
		glTranslatef(origin.X(), origin.Y(), origin.Z());
		vcg::glutSolidSphere(r, m, n); 
		glPopMatrix();

		//glEnd();
		glFlush();                                                                                                             
		glDisable(GL_LINE_SMOOTH);
	}

	void drawStickMapOnface(vcg::Point3f &origin,vcg::Point3f &dst,vcg::Point3f &normal,Color &color)
	{//////////////////////////////////////////
		// ����ͷ��
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
		vcg::glutSolidSphere(r, m, n); 
		glPopMatrix();

		glEnd();
		glFlush();                                                                                                             
		glDisable(GL_LINE_SMOOTH);
	}


	void drawArrowOnFace(vcg::Point3f &origin,vcg::Point3f &dst,vcg::Point3f &normal,Color &color)
	{//////////////////////////////////////////
		vcg::Point3f z(0,0,1);
		vcg::Point3f x_axis(1,0,0);
		vcg::Point3f y(0,1,0);	
		vcg::Point3f zero(0,0,0);

		// ��xoy���ϻ���ͷ�������Ϊdst��
		double rad = 0.15;	
		vcg::Point3f g = dst - origin;
		float wing = g.Norm()*0.3;


		// ������ת��
		vcg::Point3f axis_n_z = z ^ normal;
		vcg::Point3f axis_g_x = z;
		// ������ת�Ƕ�
		double angle_n_z = rad2angle(acos(normal*z/normal.Norm()));// n��z�ļн�	
		double angle_g_x = rad2angle(acos((g^normal)*y/g.Norm()));// g��x�ļн�

		/// ������ͷ��ֱ��
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

		// ����ͷ��
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

	// ���ֲ��������/С�Ķ��㣬������
	void showSingluarity(CFaceO& f,std::set<int>& v_Chected){
		for(int i = 0; i < 3; i++){
			CVertexO* v = f.V(i);
			if(v_Chected.count(v->Index())){// ��鶥���Ƿ��Ѿ������
				continue;
			}
			v_Chected.insert(v->Index());
			gdut_base::VertexType vt = gdut_base::countCenter(v);
			switch(vt){
			case gdut_base::SOURCE:
				drawStick(v->P() ,v->P(),gdut_base::Red);
				break;
			case gdut_base::SINK:
				drawStick(v->P() ,v->P(),gdut_base::Green);
				break;
			}
		}
	}

}

#endif