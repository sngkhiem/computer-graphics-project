#include "filter_normal.h"
#include <libqhull_r/merge_r.h>

using namespace std;
using namespace vcg;

QhullPlugin::QhullPlugin()
{
	typeList = {FP_TEST};

	for (ActionIDType tt : types())
		actionList.push_back(new QAction(filterName(tt), this));
}

QhullPlugin::~QhullPlugin()
{
}

QString QhullPlugin::pluginName() const
{
	return "FilterNormal";
}

QString QhullPlugin::filterName(ActionIDType f) const
{
	switch (f) {
	case FP_TEST:
		return QString("*** Filter Normal");
	default: assert(0); return QString();
	}
}

QString QhullPlugin::pythonFilterName(ActionIDType filterId) const
{
	switch (filterId) {
	case FP_TEST:
		return QString("Plugin name here (Python MeshLab)");
	default: assert(0); return QString();
	}
}

QString QhullPlugin::filterInfo(ActionIDType filterId) const
{
	switch (filterId) {
	case FP_TEST:
		return QString("Insert normal vector for each vertices");
	default: assert(0);
	}
	return QString("Error: Unknown Filter");
}

QhullPlugin::FilterClass QhullPlugin::getClass(const QAction* a) const
{
	switch (ID(a)) {
	case FP_TEST: return FilterClass(FilterPlugin::Remeshing);
	default: assert(0);
	}
	return FilterClass(0);
}

RichParameterList QhullPlugin::initParameterList(const QAction* action, const MeshModel& m)
{
	RichParameterList parlst;
	switch (ID(action)) {
	case FP_TEST:
		parlst.addParam(RichBool(
			"isTerrain",
			false,
			"Click if the mesh is a terrain",
			""));
	default: break; // do not add any parameter for the other filters
	}
	return parlst;
}

float calculateAngle(Point3f O, Point3f A, Point3f B)
{
	float OA = (O - A).Norm();
	float OB = (O - B).Norm();

	return acos((O - A).dot(O - B) / (2.0f * OA * OB));
}

std::map<std::string, QVariant> QhullPlugin::applyFilter(
	const QAction*           filter,
	const RichParameterList& par,
	MeshDocument&            md,
	unsigned int& /*postConditionMask*/,
	vcg::CallBackPos* /* cb*/)
{
	switch (ID(filter)) {
	case FP_TEST: {
		MeshModel& m         = *md.mm();
		//MeshModel& nm        = *md.addNewMesh("", "Mesh create by moving vertex allong normal");
		int        dim       = 3;
		int        numpoints = m.cm.vn;
		bool       isTerrain = par.getBool("isTerrain");

		m.updateDataMask(MeshModel::MM_FACEFACETOPO);
		m.updateDataMask(MeshModel::MM_VERTFACETOPO);
		m.cm.face.EnableFFAdjacency();
		m.cm.face.EnableColor();

		// reset normal of each vertex and face

		log("vertex normal %f, %f, %f",
			m.cm.face[0].V(0)->N().X(),
			m.cm.face[0].V(0)->N().Y(),
			m.cm.face[0].V(0)->N().Z());

		for (int i = 0; i < numpoints; i++) {
			m.cm.vert[i].N() = Point3f(0, 0, 0);
		}

		log("-------");

		for (int i = 0; i < m.cm.face.size(); i++) {
			m.cm.face[i].N() = Point3f(0, 0, 0);
		}

		vector<int> queue; // use to traverse though all face
		vector<int> visited(m.cm.face.size(), 0); // mark is a face is visited or not
		vector<int> ccw(m.cm.face.size(), 0); // check the vertex order store in each face 

		int         cnt = 1;
		queue.push_back(0); // push the sead face first
		visited[queue[0]] = 1; 
		ccw[queue[0]]     = 1; // assume the seed face has ccw order

		Point3f v0 = m.cm.face[queue[0]].V(0)->P();
		Point3f v1 = m.cm.face[queue[0]].V(1)->P();
		Point3f v2 = m.cm.face[queue[0]].V(2)->P();

		// (v1 - v0) x (v2 - v0)
		
		// compute normal

		float v0Angle = calculateAngle(v0, v1, v2);
		float v1Angle = calculateAngle(v1, v2, v0);
		float v2Angle = calculateAngle(v2, v0, v1);

		m.cm.face[queue[0]].V(0)->N() +=
			((v1 - v0) ^ (v2 - v0)).Scale(Point3f(v0Angle, v0Angle, v0Angle));
		m.cm.face[queue[0]].V(1)->N() +=
			((v1 - v0) ^ (v2 - v0)).Scale(Point3f(v1Angle, v1Angle, v1Angle));
		m.cm.face[queue[0]].V(2)->N() +=
			((v1 - v0) ^ (v2 - v0)).Scale(Point3f(v2Angle, v2Angle, v2Angle));

		m.cm.face[queue[0]].N() += (v1 - v0) ^ (v2 - v0);

		log("normal x: %f, y: %f, z: %f",
			((v1 - v0) ^ (v2 - v0)).X(),
			((v1 - v0) ^ (v2 - v0)).Y(),
			((v1 - v0) ^ (v2 - v0)).Z());


		for (int i = 0; i < cnt; i++) {
			CMeshO::FacePointer fp   = &(m.cm.face[queue[i]]);

			for (int ii = 0; ii < 3; ii++) {
				// find the adjacent face that does not contain vertev V(iii) of current face
				CMeshO::FacePointer fadj = fp->FFp((ii + 1) % 3);

				if (visited[fadj->Index()] == 0) {
					visited[fadj->Index()] = 1;
					queue.push_back(fadj->Index()); // add this adjacent face to
					                                //this queue if have not visited yet
					cnt++;
				}
				else {
					continue;
				}

				int index = 0; // index of the vertex of adj face that the
				               //same as V((ii + 1) % 3) of current face

				while (true) {
					if (fadj->V(index)->Index() == fp->V((ii + 1) % 3)->Index()) {
						break;
					}
					if (index < 3) {
						index++;
					}else break;
				}

				Point3f vadj0 = fadj->V(0)->P();
				Point3f vadj1 = fadj->V(1)->P();
				Point3f vadj2 = fadj->V(2)->P();

				// if fadj have the same order as current face, the vertex next
				// to C of current face (A) is different to the vertex next to C
				// of adj face (D)
				//          A * * * D
				//         * * adj *
				//        *   *   *
				//		 * cur * *
				//      B * * * C

				if (fadj->V((index + 1) % 3)->Index() != fp->V((ii + 2) % 3)->Index())
				{
					ccw[fadj->Index()] = ccw[queue[i]];
				}
				else
				{
					ccw[fadj->Index()] = 1 - ccw[queue[i]];
				}

				float   vadj0Angle   = calculateAngle(vadj0, vadj1, vadj2);
				float	vadj1Angle	 = calculateAngle(vadj1, vadj2, vadj0);
				float   vadj2Angle   = calculateAngle(vadj2, vadj0, vadj1);

				// (v1 - v0) x (v2 - v0)
				
				Point3f vertexNormal = (((vadj1 - vadj0) ^ (vadj2 - vadj0))).Normalize();

				// calculate the normal using formula N_j = Sum(a_iN_i);
				// where a_i is the angle of vertex j in face i
				// N_i is the normal of face i

				if (ccw[fadj->Index()] == 1) {
					fadj->V(0)->N() +=
						vertexNormal.Scale(Point3f(vadj0Angle, vadj0Angle, vadj0Angle));
					fadj->V(1)->N() +=
						vertexNormal.Scale(Point3f(vadj1Angle, vadj1Angle, vadj1Angle));
					fadj->V(2)->N() +=
						vertexNormal.Scale(Point3f(vadj2Angle, vadj2Angle, vadj2Angle));

					fadj->N() += (vadj1 - vadj0) ^ (vadj2 - vadj0);
				}
				else {
					fadj->V(0)->N() -=
						vertexNormal.Scale(Point3f(vadj0Angle, vadj0Angle, vadj0Angle));
					fadj->V(1)->N() -=
						vertexNormal.Scale(Point3f(vadj1Angle, vadj1Angle, vadj1Angle));
					fadj->V(2)->N() -=
						vertexNormal.Scale(Point3f(vadj2Angle, vadj2Angle, vadj2Angle));

					fadj->N() -= (vadj1 - vadj0) ^ (vadj2 - vadj0);
				}	
			}
		}
		queue.clear();
		visited.clear();
		ccw.clear();

		float baseSurfaceArea  = 0.0f;
		float draftSurfaceArea = 0.0f;
		Point3f scale            = {
            1.0f / sqrtf((float) numpoints),
            1.0f / sqrtf((float) numpoints),
            1.0f / sqrtf((float) numpoints)};

		//Normalize and make sure the normal is pointing at the righ direction

		for (int i = 0; i < numpoints; i++) {
			m.cm.vert[i].N().Normalize();
		}

		if (!isTerrain) {
			for (int i = 0; i < m.cm.face.size(); i++) {
				m.cm.face[i].N().Normalize();

				baseSurfaceArea += (m.cm.face[i].V(0)->P() - m.cm.face[i].V(1)->P())
									   .dot(m.cm.face[i].V(0)->P() - m.cm.face[i].V(2)->P());

				draftSurfaceArea +=
					(m.cm.face[i].V(0)->P() + m.cm.face[i].V(0)->N().Scale(scale) -
					 m.cm.face[i].V(1)->P() - m.cm.face[i].V(1)->N().Scale(scale))
						.dot(
							m.cm.face[i].V(0)->P() + m.cm.face[i].V(0)->N().Scale(scale) -
							m.cm.face[i].V(2)->P() - m.cm.face[i].V(2)->N().Scale(scale));
			}

			//int numOfFace = m.cm.face.size();
			//for (int i = 0; i < numOfFace; i++) {
			//	tri::Allocator<CMeshO>::AddFace(
			//		m.cm,
			//		m.cm.face[i].V(0)->P() + m.cm.face[i].V(0)->N(),
			//		m.cm.face[i].V(1)->P() + m.cm.face[i].V(1)->N(),
			//		m.cm.face[i].V(2)->P() + m.cm.face[i].V(2)->N());
			//	m.cm.face[i].C() = Color4b::LightGray;
			//	m.cm.face[numOfFace + i].C() = Color4b(0, 0, 0, 0);
			//}

			log("%f", baseSurfaceArea);
			log("%f", draftSurfaceArea);

			// if the normal of the mesh is point outward, when we move all the vertex
			// along the normal direction, the surface area must increase

			if (baseSurfaceArea > draftSurfaceArea) {
				for (int i = 0; i < numpoints; i++) {
					m.cm.vert[i].N() = -m.cm.vert[i].N();
				}

				for (int i = 0; i < m.cm.face.size(); i++) {
					m.cm.face[i].N() = -m.cm.face[i].N();
				}
			}
		}
		else {
			for (int i = 0; i < m.cm.face.size(); i++) {
				m.cm.face[i].N().Normalize();
			}
		}
		tri::Clean<CMeshO>::RemoveDuplicateVertex(m.cm);
	} break;
	default: wrongActionCalled(filter);
	}
	return std::map<std::string, QVariant>();
}
MESHLAB_PLUGIN_NAME_EXPORTER(QhullPlugin)
#include "filter_normal.h"
#include <libqhull_r/merge_r.h>

using namespace std;
using namespace vcg;

QhullPlugin::QhullPlugin()
{
	typeList = {FP_TEST};

	for (ActionIDType tt : types())
		actionList.push_back(new QAction(filterName(tt), this));
}

QhullPlugin::~QhullPlugin()
{
}

QString QhullPlugin::pluginName() const
{
	return "FilterNormal";
}

QString QhullPlugin::filterName(ActionIDType f) const
{
	switch (f) {
	case FP_TEST:
		return QString("*** Filter Normal");
	default: assert(0); return QString();
	}
}

QString QhullPlugin::pythonFilterName(ActionIDType filterId) const
{
	switch (filterId) {
	case FP_TEST:
		return QString("Plugin name here (Python MeshLab)");
	default: assert(0); return QString();
	}
}

QString QhullPlugin::filterInfo(ActionIDType filterId) const
{
	switch (filterId) {
	case FP_TEST:
		return QString("Insert normal vector for each vertices");
	default: assert(0);
	}
	return QString("Error: Unknown Filter");
}

QhullPlugin::FilterClass QhullPlugin::getClass(const QAction* a) const
{
	switch (ID(a)) {
	case FP_TEST: return FilterClass(FilterPlugin::Remeshing);
	default: assert(0);
	}
	return FilterClass(0);
}

RichParameterList QhullPlugin::initParameterList(const QAction* action, const MeshModel& m)
{
	RichParameterList parlst;
	switch (ID(action)) {
	case FP_TEST:
		parlst.addParam(RichBool(
			"isTerrain",
			false,
			"Click if the mesh is a terrain",
			""));
	default: break; // do not add any parameter for the other filters
	}
	return parlst;
}

float calculateAngle(Point3f O, Point3f A, Point3f B)
{
	float OA = (O - A).Norm();
	float OB = (O - B).Norm();

	return acos((O - A).dot(O - B) / (2.0f * OA * OB));
}

std::map<std::string, QVariant> QhullPlugin::applyFilter(
	const QAction*           filter,
	const RichParameterList& par,
	MeshDocument&            md,
	unsigned int& /*postConditionMask*/,
	vcg::CallBackPos* /* cb*/)
{
	switch (ID(filter)) {
	case FP_TEST: {
		MeshModel& m         = *md.mm();
		//MeshModel& nm        = *md.addNewMesh("", "Mesh create by moving vertex allong normal");
		int        dim       = 3;
		int        numpoints = m.cm.vn;
		bool       isTerrain = par.getBool("isTerrain");

		m.updateDataMask(MeshModel::MM_FACEFACETOPO);
		m.updateDataMask(MeshModel::MM_VERTFACETOPO);
		m.cm.face.EnableFFAdjacency();
		m.cm.face.EnableColor();

		// reset normal of each vertex and face

		log("vertex normal %f, %f, %f",
			m.cm.face[0].V(0)->N().X(),
			m.cm.face[0].V(0)->N().Y(),
			m.cm.face[0].V(0)->N().Z());

		for (int i = 0; i < numpoints; i++) {
			m.cm.vert[i].N() = Point3f(0, 0, 0);
		}

		log("-------");

		for (int i = 0; i < m.cm.face.size(); i++) {
			m.cm.face[i].N() = Point3f(0, 0, 0);
		}

		vector<int> queue; // use to traverse though all face
		vector<int> visited(m.cm.face.size(), 0); // mark is a face is visited or not
		vector<int> ccw(m.cm.face.size(), 0); // check the vertex order store in each face 

		int         cnt = 1;
		queue.push_back(0); // push the sead face first
		visited[queue[0]] = 1; 
		ccw[queue[0]]     = 1; // assume the seed face has ccw order

		Point3f v0 = m.cm.face[queue[0]].V(0)->P();
		Point3f v1 = m.cm.face[queue[0]].V(1)->P();
		Point3f v2 = m.cm.face[queue[0]].V(2)->P();

		// (v1 - v0) x (v2 - v0)
		
		// compute normal

		float v0Angle = calculateAngle(v0, v1, v2);
		float v1Angle = calculateAngle(v1, v2, v0);
		float v2Angle = calculateAngle(v2, v0, v1);

		m.cm.face[queue[0]].V(0)->N() +=
			((v1 - v0) ^ (v2 - v0)).Scale(Point3f(v0Angle, v0Angle, v0Angle));
		m.cm.face[queue[0]].V(1)->N() +=
			((v1 - v0) ^ (v2 - v0)).Scale(Point3f(v1Angle, v1Angle, v1Angle));
		m.cm.face[queue[0]].V(2)->N() +=
			((v1 - v0) ^ (v2 - v0)).Scale(Point3f(v2Angle, v2Angle, v2Angle));

		m.cm.face[queue[0]].N() += (v1 - v0) ^ (v2 - v0);

		log("normal x: %f, y: %f, z: %f",
			((v1 - v0) ^ (v2 - v0)).X(),
			((v1 - v0) ^ (v2 - v0)).Y(),
			((v1 - v0) ^ (v2 - v0)).Z());


		for (int i = 0; i < cnt; i++) {
			CMeshO::FacePointer fp   = &(m.cm.face[queue[i]]);

			for (int ii = 0; ii < 3; ii++) {
				// find the adjacent face that does not contain vertev V(iii) of current face
				CMeshO::FacePointer fadj = fp->FFp((ii + 1) % 3);

				if (visited[fadj->Index()] == 0) {
					visited[fadj->Index()] = 1;
					queue.push_back(fadj->Index()); // add this adjacent face to
					                                //this queue if have not visited yet
					cnt++;
				}
				else {
					continue;
				}

				int index = 0; // index of the vertex of adj face that the
				               //same as V((ii + 1) % 3) of current face

				while (true) {
					if (fadj->V(index)->Index() == fp->V((ii + 1) % 3)->Index()) {
						break;
					}
					if (index < 3) {
						index++;
					}else break;
				}

				Point3f vadj0 = fadj->V(0)->P();
				Point3f vadj1 = fadj->V(1)->P();
				Point3f vadj2 = fadj->V(2)->P();

				// if fadj have the same order as current face, the vertex next
				// to C of current face (A) is different to the vertex next to C
				// of adj face (D)
				//          A * * * D
				//         * * adj *
				//        *   *   *
				//		 * cur * *
				//      B * * * C

				if (fadj->V((index + 1) % 3)->Index() != fp->V((ii + 2) % 3)->Index())
				{
					ccw[fadj->Index()] = ccw[queue[i]];
				}
				else
				{
					ccw[fadj->Index()] = 1 - ccw[queue[i]];
				}

				float   vadj0Angle   = calculateAngle(vadj0, vadj1, vadj2);
				float	vadj1Angle	 = calculateAngle(vadj1, vadj2, vadj0);
				float   vadj2Angle   = calculateAngle(vadj2, vadj0, vadj1);

				// (v1 - v0) x (v2 - v0)
				
				Point3f vertexNormal = (((vadj1 - vadj0) ^ (vadj2 - vadj0))).Normalize();

				// calculate the normal using formula N_j = Sum(a_iN_i);
				// where a_i is the angle of vertex j in face i
				// N_i is the normal of face i

				if (ccw[fadj->Index()] == 1) {
					fadj->V(0)->N() +=
						vertexNormal.Scale(Point3f(vadj0Angle, vadj0Angle, vadj0Angle));
					fadj->V(1)->N() +=
						vertexNormal.Scale(Point3f(vadj1Angle, vadj1Angle, vadj1Angle));
					fadj->V(2)->N() +=
						vertexNormal.Scale(Point3f(vadj2Angle, vadj2Angle, vadj2Angle));

					fadj->N() += (vadj1 - vadj0) ^ (vadj2 - vadj0);
				}
				else {
					fadj->V(0)->N() -=
						vertexNormal.Scale(Point3f(vadj0Angle, vadj0Angle, vadj0Angle));
					fadj->V(1)->N() -=
						vertexNormal.Scale(Point3f(vadj1Angle, vadj1Angle, vadj1Angle));
					fadj->V(2)->N() -=
						vertexNormal.Scale(Point3f(vadj2Angle, vadj2Angle, vadj2Angle));

					fadj->N() -= (vadj1 - vadj0) ^ (vadj2 - vadj0);
				}	
			}
		}
		queue.clear();
		visited.clear();
		ccw.clear();

		float baseSurfaceArea  = 0.0f;
		float draftSurfaceArea = 0.0f;
		Point3f scale            = {
            1.0f / sqrtf((float) numpoints),
            1.0f / sqrtf((float) numpoints),
            1.0f / sqrtf((float) numpoints)};

		//Normalize and make sure the normal is pointing at the righ direction

		for (int i = 0; i < numpoints; i++) {
			m.cm.vert[i].N().Normalize();
		}

		if (!isTerrain) {
			for (int i = 0; i < m.cm.face.size(); i++) {
				m.cm.face[i].N().Normalize();

				baseSurfaceArea += (m.cm.face[i].V(0)->P() - m.cm.face[i].V(1)->P())
									   .dot(m.cm.face[i].V(0)->P() - m.cm.face[i].V(2)->P());

				draftSurfaceArea +=
					(m.cm.face[i].V(0)->P() + m.cm.face[i].V(0)->N().Scale(scale) -
					 m.cm.face[i].V(1)->P() - m.cm.face[i].V(1)->N().Scale(scale))
						.dot(
							m.cm.face[i].V(0)->P() + m.cm.face[i].V(0)->N().Scale(scale) -
							m.cm.face[i].V(2)->P() - m.cm.face[i].V(2)->N().Scale(scale));
			}

			//int numOfFace = m.cm.face.size();
			//for (int i = 0; i < numOfFace; i++) {
			//	tri::Allocator<CMeshO>::AddFace(
			//		m.cm,
			//		m.cm.face[i].V(0)->P() + m.cm.face[i].V(0)->N(),
			//		m.cm.face[i].V(1)->P() + m.cm.face[i].V(1)->N(),
			//		m.cm.face[i].V(2)->P() + m.cm.face[i].V(2)->N());
			//	m.cm.face[i].C() = Color4b::LightGray;
			//	m.cm.face[numOfFace + i].C() = Color4b(0, 0, 0, 0);
			//}

			log("%f", baseSurfaceArea);
			log("%f", draftSurfaceArea);

			// if the normal of the mesh is point outward, when we move all the vertex
			// along the normal direction, the surface area must increase

			if (baseSurfaceArea > draftSurfaceArea) {
				for (int i = 0; i < numpoints; i++) {
					m.cm.vert[i].N() = -m.cm.vert[i].N();
				}

				for (int i = 0; i < m.cm.face.size(); i++) {
					m.cm.face[i].N() = -m.cm.face[i].N();
				}
			}
		}
		else {
			for (int i = 0; i < m.cm.face.size(); i++) {
				m.cm.face[i].N().Normalize();
			}
		}
		tri::Clean<CMeshO>::RemoveDuplicateVertex(m.cm);
	} break;
	default: wrongActionCalled(filter);
	}
	return std::map<std::string, QVariant>();
}
MESHLAB_PLUGIN_NAME_EXPORTER(QhullPlugin)
