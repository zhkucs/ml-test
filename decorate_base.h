/****************************************************************************
* MeshLab                                                           o o     *
* A versatile mesh processing toolbox                             o     o   *
*                                                                _   O  _   *
* Copyright(C) 2005                                                \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/

#ifndef EXTRADECORATEPLUGIN_H
#define EXTRADECORATEPLUGIN_H

#include <common/interfaces.h>
#include <wrap/gui/coordinateframe.h>
#include "colorhistogram.h"
#include "alg_base.h"



class QGLShaderProgram;
typedef vcg::ColorHistogram<float>  CHist;

typedef std::pair<vcg::Point3f,vcg::Color4b> PointPC; // this type is used to have a simple coord+color pair to rapidly draw non manifold faces

class ExtraMeshDecoratePlugin : public QObject, public MeshDecorateInterface
{
  Q_OBJECT
    MESHLAB_PLUGIN_IID_EXPORTER(MESH_DECORATE_INTERFACE_IID)
  Q_INTERFACES(MeshDecorateInterface)
  QString decorationName(FilterIDType filter) const;
  QString decorationInfo(FilterIDType filter) const;

  enum {
    DP_SHOW_NORMALS,
    DP_SHOW_VERT,
    DP_SHOW_EDGE,
    DP_SHOW_NON_FAUX_EDGE,
    DP_SHOW_BOUNDARY,
    DP_SHOW_NON_MANIF_EDGE,
    DP_SHOW_NON_MANIF_VERT,
    DP_SHOW_BOX_CORNERS,
    DP_SHOW_AXIS,
    DP_SHOW_QUOTED_BOX,
    DP_SHOW_LABEL,
    DP_SHOW_QUALITY_HISTOGRAM,
    DP_SHOW_QUALITY_CONTOUR,
    DP_SHOW_CAMERA,
    DP_SHOW_TEXPARAM,
    DP_SHOW_BOUNDARY_TEX,
    DP_SHOW_SELECTED_MESH,
	DP_SHOW_GRADIENT,
	DP_MEAN_CURVATURE,
	DP_HHD_CURLFREE,
	DP_HHD_DIVFREE,
	DP_HHD_HARMONEY
	
	/*,
    DP_SHOW_SELECTED_FACE,
    DP_SHOW_SELECTED_VERT*/
  };


private:
	float niceRound2(float value,float base);
	float niceRound(float value);

  void	drawQuotedLine(const vcg::Point3d &a,const vcg::Point3d &b,float aVal, float bVal,float tickDist,QPainter *painter, QFont qf,float angle =0,bool rightAlign=false);

	void	chooseX(vcg::Box3f &box,double *modelview,double *projection,GLint *viewport,vcg::Point3d &x1,vcg::Point3d &x2);
	void	chooseY(vcg::Box3f &box,double *modelview,double *projection,GLint *viewport,vcg::Point3d &y1,vcg::Point3d &y2);
	void	chooseZ(vcg::Box3f &box,double *modelview,double *projection,GLint *viewport,vcg::Point3d &z1,vcg::Point3d &z2);
  void drawHistogram(QGLWidget *gla, CHist &ch);

  bool countOnce;
  std::map<int,vcg::Point3f> m_kexi;
  Eigen::VectorXd m_x;
  std::map<int,vcg::Point3f> m_df;

  std::set<vcg::Point3f> s_source;
  std::set<vcg::Point3f> s_sink;

  std::map<std::string,ShowParameter*> paraMap;
public:

  ExtraMeshDecoratePlugin()
  {
    typeList <<
                DP_SHOW_VERT <<
                //DP_SHOW_EDGE <<
                DP_SHOW_NON_FAUX_EDGE <<
                DP_SHOW_BOUNDARY <<
                DP_SHOW_NON_MANIF_EDGE <<
                DP_SHOW_NON_MANIF_VERT <<
                DP_SHOW_NORMALS <<
                DP_SHOW_QUALITY_HISTOGRAM <<
                DP_SHOW_QUALITY_CONTOUR <<
                DP_SHOW_BOX_CORNERS <<
                DP_SHOW_AXIS <<
                DP_SHOW_QUOTED_BOX <<
                DP_SHOW_LABEL <<
                DP_SHOW_CAMERA <<
                DP_SHOW_TEXPARAM <<
                DP_SHOW_SELECTED_MESH <<
                /*DP_SHOW_SELECTED_FACE <<
                DP_SHOW_SELECTED_VERT <<*/
                DP_SHOW_BOUNDARY_TEX<<
				// MY MENU
				DP_SHOW_GRADIENT<<
				DP_MEAN_CURVATURE<<
				DP_HHD_CURLFREE<<
				DP_HHD_DIVFREE<<
				DP_HHD_HARMONEY
				;

    FilterIDType tt;
    foreach(tt , types())
    {
      actionList << new QAction(decorationName(tt), this);
      /*if(tt==DP_SHOW_SELECTED_VERT)
        actionList.last()->setIcon(QIcon(":/images/selected_vert.png"));
      if(tt==DP_SHOW_SELECTED_FACE)
        actionList.last()->setIcon(QIcon(":/images/selected_face.png"));*/
    }

  QAction *ap;
  foreach(ap,actionList){
    ap->setCheckable(true);
  }

  countOnce = false;

  initMap();

}
  void initMap();
  void DrawBBoxCorner(MeshModel &m, bool absBBoxFlag=true);
  void DrawQuotedBox(MeshModel &m,QPainter *gla, QFont qf);
  void DrawVertLabel(MeshModel &m, QPainter *gla);
  void DrawEdgeLabel(MeshModel &m, QPainter *gla);
  void DrawFaceLabel(MeshModel &m, QPainter *gla);
  void DisplayCamera(MeshModel *m, vcg::Shotf &ls, int cameraSourceId, QPainter *painter, QFont qf);
  void DrawCamera(MeshModel *m, vcg::Shotf &ls, vcg::Color4b camcolor, vcg::Matrix44f &currtr, RichParameterSet *rm, QPainter *painter, QFont qf);
  void PlaceTexParam(int TexInd, int TexNum);
  void DrawTexParam(MeshModel &m, GLArea *gla, QPainter *painter, RichParameterSet *, QFont qf);
  void DrawColorHistogram(CHist &ch, GLArea *gla, QPainter *painter, RichParameterSet *, QFont qf);
  void DrawLineVector(std::vector<PointPC> &EV);
  void DrawTriVector(std::vector<PointPC> &EV);
  void DrawDotVector(std::vector<PointPC> &EV, float basesize=4.0);



  void decorateDoc(QAction *a, MeshDocument &md, RichParameterSet *, GLArea *gla, QPainter *painter, GLLogStream &_log);
  void decorateMesh(QAction *a, MeshModel &md, RichParameterSet *, GLArea *gla, QPainter *painter, GLLogStream &_log);
  bool startDecorate(QAction * /*mode*/, MeshModel &/*m*/, RichParameterSet *, GLArea * /*parent*/);
  void endDecorate(QAction * /*mode*/, MeshModel &/*m*/, RichParameterSet *, GLArea * /*parent*/);
  bool startDecorate(QAction * /*mode*/, MeshDocument &/*m*/, RichParameterSet *, GLArea * /*parent*/);
  bool isDecorationApplicable(QAction *action, const MeshModel& m, QString &ErrorMessage) const;
  int getDecorationClass(QAction */*action*/) const;
  void initGlobalParameterSet(QAction *, RichParameterSet &/*globalparam*/);

  inline QString CameraScaleParam() const    { return  "MeshLab::Decoration::CameraRenderScaleType" ; }
  inline QString FixedScaleParam() const     { return  "MeshLab::Decoration::CameraFixedScaleParam" ; }

  inline QString ShowCameraDetails() const     { return  "MeshLab::Decoration::CameraShowCameraDetails" ; }

  inline QString ShowMeshCameras() const     { return  "MeshLab::Decoration::ShowMeshCameras" ; }
  inline QString ShowRasterCameras() const   { return  "MeshLab::Decoration::ShowRasterCameras" ; }

  inline QString ShowNonRegular() const     { return  "MeshLab::Decoration::ShowNonRegular" ; }
  inline QString ShowSeparatrix() const     { return  "MeshLab::Decoration::ShowSeparatrix" ; }
  inline QString ShowBorderFlag() const     { return  "MeshLab::Decoration::ShowBorderFlag" ; }

  inline QString TextureStyleParam() const   { return  "MeshLab::Decoration::TextureStyle" ; }
  inline QString TextureFaceColorParam() const   { return  "MeshLab::Decoration::TextureFaceColor" ; }
  inline QString VertDotSizeParam() const    { return  "MeshLab::Decoration::VertDotSize" ; }

  inline QString HistBinNumParam() const     { return  "MeshLab::Decoration::HistBinNumParam" ; }
  inline QString HistFixedParam() const   { return  "MeshLab::Decoration::UseFixedHistParam" ; }
  inline QString HistFixedMinParam() const   { return  "MeshLab::Decoration::FixedHistMinParam" ; }
  inline QString HistFixedMaxParam() const   { return  "MeshLab::Decoration::FixedHistMaxParam" ; }
  inline QString HistFixedWidthParam() const { return  "MeshLab::Decoration::FixedHistWidthParam" ; }
  inline QString HistAreaParam() const { return  "MeshLab::Decoration::AreaHistParam" ; }
  inline QString HistTypeParam() const { return  "MeshLab::Decoration::HistType" ; }

  inline QString ShowContourFreq() const { return  "MeshLab::Decoration::ShowContourFreq" ; }
  inline QString ShowContourAlpha() const { return  "MeshLab::Decoration::ShowContourAlpha" ; }
  inline QString ShowContourWidth() const { return  "MeshLab::Decoration::ShowContourWidth" ; }
  inline QString ShowContourRamp() const { return  "MeshLab::Decoration::ShowContourRamp" ; }

  inline QString LabelVertFlag() const { return  "MeshLab::Decoration::LabelVertFlag" ; }
  inline QString LabelEdgeFlag() const { return  "MeshLab::Decoration::LabelEdgeFlag" ; }
  inline QString LabelFaceFlag() const { return  "MeshLab::Decoration::LabelFaceFlag" ; }

  inline QString NormalLength() const { return  "MeshLab::Decoration::NormalLength" ; }
  inline QString NormalVertFlag() const { return  "MeshLab::Decoration::NormalVertFlag" ; }
  inline QString NormalFaceFlag() const { return  "MeshLab::Decoration::NormalFaceFlag" ; }
  inline QString NormalCurvFlag() const { return  "MeshLab::Decoration::NormalCurvFlag" ; }

  inline QString BBAbsParam() const { return  "MeshLab::Decoration::BBAbs" ; }

  inline vcg::Point3f standardize(vcg::Point3f& start,vcg::Point3f& end,double max)
{
	vcg::Point3f diff = (start - end);	
	double r = max / diff.Norm();
	return start + diff* r;
}

  //void drawArrow(vcg::Point3f& start,vcg::Point3f& end,vcg::Point3f &normal,float*color);

signals:
  void askViewerShot(QString);

public slots:
  void  setValue(QString name, vcg::Shotf val);

private:
  vcg::Shotf curShot;

  QMap<MeshModel *, QGLShaderProgram *> contourShaderProgramMap;
};



#endif
