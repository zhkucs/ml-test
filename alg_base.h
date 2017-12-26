#ifndef ALG_BASE_H
#define ALG_BASE_H

#include <eigen/Sparse>
#include <common/meshmodel.h>
#include <common/interfaces.h>
#include <vector>
#include <math.h>
using namespace std;
using namespace vcg;
//typedef Eigen::Triplet<double> T;
//typedef Eigen::SparseMatrix<double> SparseMatrixType;
#define EPSLON 0.001
#define PI 3.141592653589
#define rad2angle( r) r*180/PI;
namespace gdut_base{
	//�������� ���ڶ�pair����hash
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
		
	}Red(1,0,0),Green(0,1,0),Blue(0,0,1);

	void countNablaOfFace(CFaceO&f,double s0,double s1,double s2,vcg::Point3f&result){
				vcg::Point3f p_i = f.P0(0); 
				vcg::Point3f p_j = f.P0(1); 
				vcg::Point3f p_k = f.P0(2); 

				vcg::Point3f vij=f.P0(1)-f.P0(0);  
				vcg::Point3f vik=f.P0(2)-f.P0(0);
				float area = (vij ^ vik).Norm()/2;

				vcg::Point3f normalf=(vij ^ vik).Normalize(); 

				vcg::Point3f phi_j = (p_i-p_j)^normalf/(2*area);
				vcg::Point3f phi_k = (p_j-p_i)^normalf/(2*area);

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
		float double_area = (vij ^ vik).Norm()/2;
		vcg::Point3f normalf=(vij ^ vik).Normalize(); 
		result= (p_j-p_k)^normalf/double_area;
	}

	float countArea(CFaceO& f){
		vcg::Point3f p_i = f.P(0);
		vcg::Point3f p_j = f.P(1);
		vcg::Point3f p_k = f.P(2);
		vcg::Point3f vij=p_j-p_i;  
		vcg::Point3f vik=-p_k-p_i;
		return (vij ^ vik).Norm()/2;
	}


	void drawArrow(vcg::Point3f &origin,vcg::Point3f &dst,vcg::Point3f &normal,Color &color)
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

}

#endif