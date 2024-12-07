#include "filter_idt.h"
#include <Eigen/Sparse>
#include <cmath>
#include <libqhull_r/merge_r.h>
#include <math.h>
#include <fstream>
#include <iostream>

using namespace std;
using namespace vcg;
using namespace Eigen;

const double PI = 3.141592653589793238462643383279502884L;

QhullPlugin::QhullPlugin()
{
	typeList = {FP_IDT};

	for (ActionIDType tt : types())
		actionList.push_back(new QAction(filterName(tt), this));
}

QhullPlugin::~QhullPlugin()
{
}

QString QhullPlugin::pluginName() const
{
	return "FilterIDT";
}

QString QhullPlugin::filterName(ActionIDType f) const
{
	switch (f) {
	case FP_IDT: return QString("Intrinsic Delaunay Triangulation");
	default: assert(0); return QString();
	}
}

QString QhullPlugin::pythonFilterName(ActionIDType filterId) const
{
	switch (filterId) {
	case FP_IDT: return QString("Plugin name here (Python MeshLab)");
	default: assert(0); return QString();
	}
}

QString QhullPlugin::filterInfo(ActionIDType filterId) const
{
	switch (filterId) {
	case FP_IDT: return QString("test");
	default: assert(0);
	}
	return QString("Error: Unknown Filter");
}

QhullPlugin::FilterClass QhullPlugin::getClass(const QAction* a) const
{
	switch (ID(a)) {
	case FP_IDT: return FilterClass(FilterPlugin::Remeshing);
	default: assert(0);
	}
	return FilterClass(0);
}

RichParameterList QhullPlugin::initParameterList(const QAction* action, const MeshModel& m)
{
	RichParameterList parlst;
	switch (ID(action)) {
	case FP_IDT:
		parlst.addParam(RichBool(
			"heatMethodBefore",
			false,
			"Heat method before IDT",
			"Calculuate Geodesic Distance from source vertex to all vertex using Heat method."));
		parlst.addParam(RichBool(
			"heatMethodAfter",
			false,
			"Heat method after IDT",
			"Calculuate Geodesic Distance from source vertex to all vertex using Heat method."));
		parlst.addParam(
			RichBool("showIDT", false, "Show new mesh after IDT", "Using IDT to remeshing. "));
		parlst.addParam(RichInt(
			"sourceVertex",
			(int) 0,
			"The source vertex for Heat Method",
			"Choosing source vertex for Heat Method"));
	default: break;
	}
	return parlst;
}

void constructF(CMeshO& cm, vector<vector<int>>& faces)
{
	for (int i = 0; i < cm.face.size(); i++) {
		faces[i][0] = cm.face[i].V(0)->Index();
		faces[i][1] = cm.face[i].V(1)->Index();
		faces[i][2] = cm.face[i].V(2)->Index();
	}
}

void constructL(CMeshO& cm, vector<vector<int>> faces, vector<vector<double>>& lengths)
{
	for (int i = 0; i < faces.size(); i++) {
		for (int j = 0; j < 3; j++) {
			int firstPointIdx  = faces[i][j];
			int secondPointIdx = faces[i][(j + 1) % 3];

			double dx = (double) cm.vert[secondPointIdx].P().X() - cm.vert[firstPointIdx].P().X();
			double dy = (double) cm.vert[secondPointIdx].P().Y() - cm.vert[firstPointIdx].P().Y();
			double dz = (double) cm.vert[secondPointIdx].P().Z() - cm.vert[firstPointIdx].P().Z();
			double length = sqrtl(dx * dx + dy * dy + dz * dz);
			lengths[i][j] = length;
		}
	}
}

pair<int, int> other(const vector<vector<vector<int>>>& G, const pair<int, int>& fs)
{
	return {G[fs.first][fs.second][0], G[fs.first][fs.second][1]};
}

void glue_together(vector<vector<vector<int>>>& G, pair<int, int> fs0, pair<int, int> fs1)
{
	G[fs0.first][fs0.second][0] = fs1.first;
	G[fs0.first][fs0.second][1] = fs1.second;
	G[fs1.first][fs1.second][0] = fs0.first;
	G[fs1.first][fs1.second][1] = fs0.second;
}

void constructG(vector<vector<int>>& faces, vector<vector<vector<int>>>& G)
{
	int n_faces = faces.size();

	if (n_faces == 0) {
		throw runtime_error("No faces provided. Mesh is empty.");
	}

	int                               n_sides = 3 * n_faces;
	vector<tuple<int, int, int, int>> S(n_sides);

	for (int f = 0; f < n_faces; f++) {
		for (int s = 0; s < 3; s++) {
			int i = faces[f][s];
			int j = faces[f][(s + 1) % 3];

			if (i < 0 || j < 0) {
				throw runtime_error("Invalid vertex index in face.");
			}

			S[f * 3 + s] = {min(i, j), max(i, j), f, s};
		}
	}

	sort(S.begin(), S.end());

	G.resize(n_faces, vector<vector<int>>(3, {-1, -1}));

	for (int p = 0; p < n_sides;) {
		if (p + 1 < n_sides && get<0>(S[p]) == get<0>(S[p + 1]) &&
			get<1>(S[p]) == get<1>(S[p + 1])) {
			pair<int, int> fs0 = {get<2>(S[p]), get<3>(S[p])};
			pair<int, int> fs1 = {get<2>(S[p + 1]), get<3>(S[p + 1])};

			glue_together(G, fs0, fs1);
			p += 2;
		}
		else {
			int f   = get<2>(S[p]);
			int s   = get<3>(S[p]);
			G[f][s] = {-1, -1};
			p++;
		}
	}
}

double face_area(vector<vector<double>>& lengths, int f)
{
	double l_a = lengths[f][0];
	double l_b = lengths[f][1];
	double l_c = lengths[f][2];

	double s = (l_a + l_b + l_c) / 2.0;
	double d = s * (s - l_a) * (s - l_b) * (s - l_c);
	return sqrtl(d);
}

double surface_area(vector<vector<int>>& faces, vector<vector<double>>& lengths)
{
	double areaTotal = 0.0;
	for (int f = 0; f < faces.size(); f++) {
		areaTotal += face_area(lengths, f);
	}
	return areaTotal;
}

double oppositeCornerAngle(vector<vector<double>>& lengths, pair<int, int> fs)
{
	double l_a = lengths[fs.first][fs.second];
	double l_b = lengths[fs.first][(fs.second + 1) % 3];
	double l_c = lengths[fs.first][((fs.second + 1) % 3 + 1) % 3];

	double d = (l_b * l_b + l_c * l_c - l_a * l_a) / (2.0L * l_b * l_c);

	return acosl(d);
}

double
diagonal_length(vector<vector<vector<int>>>& G, vector<vector<double>>& lengths, pair<int, int> fs)
{
	pair<int, int> fs_opp  = other(G, fs);
	double         u       = lengths[fs.first][((fs.second + 1) % 3 + 1) % 3];
	double         v       = lengths[fs_opp.first][(fs_opp.second + 1) % 3];
	double         theta_A = oppositeCornerAngle(lengths, {fs.first, (fs.second + 1) % 3});
	double         theta_B =
		oppositeCornerAngle(lengths, {fs_opp.first, ((fs_opp.second + 1) % 3 + 1) % 3});
	double cosThetaSum = cosl(theta_A + theta_B);
	double d           = sqrtl(u * u + v * v - 2.0L * u * v * cosThetaSum);
	return d;
}

pair<int, int> relabel(
	pair<int, int> s,
	pair<int, int> s2,
	pair<int, int> s3,
	pair<int, int> s4,
	pair<int, int> s5,
	int            f0,
	int            f1)
{
	if (s == s2)
		return {f1, 2};
	if (s == s3)
		return {f0, 1};
	if (s == s4)
		return {f0, 2};
	if (s == s5)
		return {f1, 1};
	return s;
}

bool isDelaunay(vector<vector<vector<int>>>& G, vector<vector<double>>& lengths, pair<int, int> fs)
{
	pair<int, int> fs_opp = other(G, fs);

	double theta_A = oppositeCornerAngle(lengths, fs);
	double theta_B = oppositeCornerAngle(lengths, fs_opp);

	double EPS = 1e-5;
	return (theta_A + theta_B) <= (acosl(-1.0L) + EPS);
}

pair<int, int> flipEdge(
	vector<vector<int>>&         faces,
	vector<vector<vector<int>>>& G,
	vector<vector<double>>&      lengths,
	pair<int, int>&              s0)
{
	pair<int, int> s1 = other(G, s0);

	pair<int, int> s2 = {s0.first, (s0.second + 1) % 3};
	pair<int, int> s3 = {s0.first, ((s0.second + 1) % 3 + 1) % 3};
	pair<int, int> s4 = {s1.first, (s1.second + 1) % 3};
	pair<int, int> s5 = {s1.first, ((s1.second + 1) % 3 + 1) % 3};

	pair<int, int> s6 = other(G, s2);
	pair<int, int> s7 = other(G, s3);
	pair<int, int> s8 = other(G, s4);
	pair<int, int> s9 = other(G, s5);

	int v0 = faces[s0.first][s0.second];
	int v1 = faces[s2.first][s2.second];
	int v2 = faces[s3.first][s3.second];
	int v3 = faces[s5.first][s5.second];

	int f0 = s0.first;
	int f1 = s1.first;

	double l2 = lengths[s2.first][s2.second];
	double l3 = lengths[s3.first][s3.second];
	double l4 = lengths[s4.first][s4.second];
	double l5 = lengths[s5.first][s5.second];

	double newLength = diagonal_length(G, lengths, s0);

	faces[f0][0] = v3;
	faces[f0][1] = v2;
	faces[f0][2] = v0;
	faces[f1][0] = v2;
	faces[f1][1] = v3;
	faces[f1][2] = v1;

	s6 = relabel(s6, s2, s3, s4, s5, f0, f1);
	s7 = relabel(s7, s2, s3, s4, s5, f0, f1);
	s8 = relabel(s8, s2, s3, s4, s5, f0, f1);
	s9 = relabel(s9, s2, s3, s4, s5, f0, f1);

	glue_together(G, {f0, 0}, {f1, 0});
	glue_together(G, {f0, 1}, s7);
	glue_together(G, {f0, 2}, s8);
	glue_together(G, {f1, 1}, s9);
	glue_together(G, {f1, 2}, s6);

	lengths[f0][0] = newLength;
	lengths[f0][1] = l3;
	lengths[f0][2] = l4;
	lengths[f1][0] = newLength;
	lengths[f1][1] = l5;
	lengths[f1][2] = l2;

	return {f0, 0};
}

int flipToDelaunay(
	vector<vector<int>>&         faces,
	vector<vector<vector<int>>>& G,
	vector<vector<double>>&      lengths
)
{
	int cnt = 0;
	vector<bool>          refinedFace(faces.size(), false);
	vector<vector<bool>>  vis(faces.size(), vector<bool>(faces[0].size(), 0));
	queue<pair<int, int>> process;
	for (int f = 0; f < faces.size(); f++) {
		for (int s = 0; s < 3; s++) {
			if (G[f][s][0] == -1 && G[f][s][1] == -1)
				continue;
			process.push({f, s});
			vis[f][s] = true;
		}
	}

	while (!process.empty()) {
		pair<int, int> fs = process.front();
		qDebug("%i, %i", fs.first, fs.second);
		process.pop();
		vis[fs.first][fs.second] = false;

		if (!isDelaunay(G, lengths, fs)) {
			cnt++;

			fs = flipEdge(faces, G, lengths, fs);


			if (!vis[fs.first][(fs.second + 1) % 3] &&
				!isDelaunay(G, lengths, {fs.first, (fs.second + 1) % 3})) {
				vis[fs.first][(fs.second + 1) % 3] = true;
				process.push({fs.first, (fs.second + 1) % 3});
			}
			if (!vis[fs.first][((fs.second + 1) % 3 + 1) % 3] &&
				!isDelaunay(G, lengths, {fs.first, ((fs.second + 1) % 3 + 1) % 3})) {
				vis[fs.first][((fs.second + 1) % 3 + 1) % 3] = true;
				process.push({fs.first, ((fs.second + 1) % 3 + 1) % 3});
			}
			if (!vis[other(G, fs).first][(other(G, fs).second + 1) % 3] &&
				!isDelaunay(G, lengths, {other(G, fs).first, (other(G, fs).second + 1) % 3})) {
				vis[other(G, fs).first][(other(G, fs).second + 1) % 3] = true;
				process.push({other(G, fs).first, (other(G, fs).second + 1) % 3});
			}
			if (!vis[other(G, fs).first][((other(G, fs).second + 1) % 3 + 1) % 3] &&
				!isDelaunay(
					G, lengths, {other(G, fs).first, ((other(G, fs).second + 1) % 3 + 1) % 3})) {
				vis[other(G, fs).first][((other(G, fs).second + 1) % 3 + 1) % 3] = true;
				process.push({other(G, fs).first, ((other(G, fs).second + 1) % 3 + 1) % 3});
			}

		}
	}
	return cnt;
}

SparseMatrix<double>
buildCotanLaplacian(CMeshO& cm, vector<vector<int>>& faces, vector<vector<double>>& lengths)
{
	int                     N = cm.vert.size();
	SparseMatrix<double>    L(N, N);
	vector<Triplet<double>> triplets;

	for (int f = 0; f < faces.size(); f++) {
		for (int s = 0; s < 3; s++) {
			int i = faces[f][s];
			int j = faces[f][(s + 1) % 3];

			double opp_theta    = oppositeCornerAngle(lengths, {f, s});
			double opp_cotan    = 1.0 / tan(opp_theta);
			double cotan_weight = 0.5 * opp_cotan;

			triplets.emplace_back(i, j, -cotan_weight);
			triplets.emplace_back(j, i, -cotan_weight);
			triplets.emplace_back(i, i, cotan_weight);
			triplets.emplace_back(j, j, cotan_weight);
		}
	}

	L.setFromTriplets(triplets.begin(), triplets.end());
	return L;
}

SparseMatrix<double>
buildLumpedMass(CMeshO& cm, vector<vector<int>>& faces, vector<vector<double>> lengths)
{
	int                     N = cm.vert.size();
	SparseMatrix<double>    M(N, N);
	vector<Triplet<double>> triplets;

	for (int f = 0; f < faces.size(); f++) {
		double area = face_area(lengths, f);
		for (int s = 0; s < 3; s++) {
			int i = faces[f][s];
			triplets.emplace_back(i, i, area / 3.0L);
		}
	}

	M.setFromTriplets(triplets.begin(), triplets.end());
	return M;
}

Vector2d edgeInFaceBasis(vector<vector<double>>& lengths, pair<int, int> fs)
{
	int f = fs.first;
	int s = fs.second;

	double theta = oppositeCornerAngle(lengths, {f, 1});

	Matrix<double, 3, 2> localVertPos;
	localVertPos << 0.0, 0.0, lengths[f][0], 0.0, cosl(theta) * lengths[f][2],
		sinl(theta) * lengths[f][2];

	Vector2d edgeVec = localVertPos.row((s + 1) % 3) - localVertPos.row(s);
	return edgeVec;
}

vector<Vector2d>
evaluateGradientAtFaces(vector<vector<int>>& faces, vector<vector<double>> lengths, VectorXd& x)
{
	int              n_faces = faces.size();
	vector<Vector2d> grads(n_faces);

	for (int f = 0; f < n_faces; f++) {
		Vector2d face_grad(0.0, 0.0);

		for (int s = 0; s < 3; s++) {
			int      i       = faces[f][s];
			Vector2d edgeVec = edgeInFaceBasis(lengths, {f, (s + 1) % 3});
			Vector2d edgeVecRot(-edgeVec.y(), edgeVec.x());

			face_grad += x[i] * edgeVecRot;
		}

		double area = face_area(lengths, f);
		face_grad /= (2.0 * area);

		grads[f] = face_grad;
	}

	return grads;
}

vector<double> evaluateDivergenceAtVertices(
	CMeshO&                 cm,
	vector<vector<int>>&    faces,
	vector<vector<double>>& lengths,
	vector<Vector2d>&       v)
{
	vector<double> divs(cm.vert.size(), 0.0);

	for (int f = 0; f < faces.size(); f++) {
		Vector2d gradVec = v[f];

		for (int s = 0; s < 3; s++) {
			int i = faces[f][s];
			int j = faces[f][(s + 1) % 3];

			Vector2d edgeVec = edgeInFaceBasis(lengths, {f, s});

			double opp_theta    = oppositeCornerAngle(lengths, {f, s});
			double opp_cotan    = 1.0L / tanl(opp_theta);
			double cotan_weight = 0.5L * opp_cotan;

			double div_contrib = cotan_weight * edgeVec.dot(gradVec);

			divs[i] += div_contrib;
			divs[j] -= div_contrib;
		}
	}

	return divs;
}

void heatMethod(
	CMeshO&                 cm,
	vector<vector<int>>&    faces,
	vector<vector<double>>& lengths,
	int                     sourceVertex,
	vector<double>&         distances)
{
	SparseMatrix<double> L = buildCotanLaplacian(cm, faces, lengths);
	SparseMatrix<double> M = buildLumpedMass(cm, faces, lengths);

	double meanEdgeLength  = 0.0;
	double totalEdgeLength = 0.0;
	int    edgeCount       = 0;

	for (const auto& faceLengths : lengths) {
		for (double len : faceLengths) {
			totalEdgeLength += len;
			++edgeCount;
		}
	}
	meanEdgeLength = (double) totalEdgeLength / edgeCount;

	double shortTime = meanEdgeLength * meanEdgeLength;

	SparseMatrix<double> H = M + shortTime * L;

	VectorXd initRHS      = VectorXd::Zero(L.rows());
	initRHS[sourceVertex] = 1.0L;

	SparseLU<SparseMatrix<double>> solver;
	solver.compute(H);
	if (solver.info() != Success) {
		throw std::runtime_error("Failed to decompose heat operator matrix.");
	}
	VectorXd heat = solver.solve(initRHS);
	if (solver.info() != Success) {
		throw std::runtime_error("Failed to solve heat equation.");
	}

	vector<Vector2d> gradients = evaluateGradientAtFaces(faces, lengths, heat);
	for (auto& grad : gradients) {
		grad.normalize();
	}

	vector<double> divergence = evaluateDivergenceAtVertices(cm, faces, lengths, gradients);

	VectorXd div = VectorXd::Zero(divergence.size());
	for (int i = 0; i < divergence.size(); i++) {
		div[i] = divergence[i];
	}
	SparseMatrix<double> regularizedL = L;
	for (int k = 0; k < L.rows(); k++) {
		regularizedL.coeffRef(k, k) += 1e-6;
	}

	BiCGSTAB<SparseMatrix<double>> altSolver;
	altSolver.compute(regularizedL);
	if (altSolver.info() != Success) {
		throw std::runtime_error("Failed to decompose regularized Laplacian matrix with BiCGSTAB.");
	}
	VectorXd dist = altSolver.solve(div);
	if (altSolver.info() != Success) {
		throw std::runtime_error("Failed to solve Poisson equation with BiCGSTAB.");
	}

	distances.resize(dist.size());
	double distSrc = dist[sourceVertex];
	for (int i = 0; i < dist.size(); i++)
		dist[i] -= distSrc;
	for (int i = 0; i < dist.size(); ++i) {
		distances[i] = (dist[i]);
	}
}

std::map<std::string, QVariant> QhullPlugin::applyFilter(
	const QAction*           filter,
	const RichParameterList& par,
	MeshDocument&            md,
	unsigned int& /*postConditionMask*/,
	vcg::CallBackPos* /* cb*/)
{
	switch (ID(filter)) {
	case FP_IDT: {
		ofstream file, file2;
		file2.open("D:/CG/log2.txt");
		file.open("D:/CG/log.txt");
		file.clear();
		file2.clear();
		MeshModel&                  m = *md.mm();
		vector<vector<int>>         faces(m.cm.face.size(), vector<int>(3, 0));
		vector<vector<double>>      lengths(m.cm.face.size(), vector<double>(3, 0.0));
		vector<vector<vector<int>>> G(faces.size(), vector<vector<int>>(3, vector<int>(2, 0)));

		bool heatMethodBefore = par.getBool("heatMethodBefore");
		bool heatMethodAfter  = par.getBool("heatMethodAfter");
		bool showIDT          = par.getBool("showIDT");
		int  sourceVertex     = par.getInt("sourceVertex");

		constructF(m.cm, faces);
		constructL(m.cm, faces, lengths);
		constructG(faces, G);

		if (heatMethodBefore) {
			vector<double> distances(m.cm.vert.size());
			heatMethod(m.cm, faces, lengths, sourceVertex, distances);
			file2 << faces.size() << "\n";
			for (int i = 0; i < faces.size(); i++) {
				file2 << faces[i][0] << " " << faces[i][1] << " " << faces[i][2] << "\n";
			}
			file2 << m.cm.vert.size() << "\n" << sourceVertex << "\n";
			for (int i = 0; i < m.cm.vert.size(); i++)
				file2 << distances[i] << "\n";
			for (int i = 0; i < m.cm.vert.size(); i++)
				file2 << m.cm.vert[i].P().X() << " " << m.cm.vert[i].P().Y() << " " << m.cm.vert[i].P().Z() << "\n";
			file2.close();
			MeshModel& nm1 = *md.addNewMesh("", "Heat Method Before IDT");
			tri::Append<CMeshO, CMeshO>::MeshCopy(nm1.cm, m.cm);
			nm1.cm.face.EnableColor();

			vector<array<double, 3>> colors(distances.size());

			double minDist = *min_element(distances.begin(), distances.end());
			double maxDist = *max_element(distances.begin(), distances.end());
			double range   = maxDist - minDist;

			double radiusStep = 0.05;

			for (size_t i = 0; i < nm1.cm.face.size(); ++i) {
				CMeshO::FaceType& face = nm1.cm.face[i];

				double avgDistance = 0.0;

				for (int j = 0; j < 3; ++j) {
					int vertexIndex = face.V(j)->Index();
					avgDistance += distances[vertexIndex];
				}
				avgDistance /= 3.0;

				double normalized = (avgDistance - minDist) / range;
				int    band       = static_cast<int>(floor(normalized / radiusStep));

				double green = (band % 2 == 0) ? 255.0 : 50.0;

				face.C() = vcg::Color4b(0, static_cast<unsigned char>(green), 0, 255);
			}

			nm1.updateDataMask(MeshModel::MM_FACECOLOR);

		}

		log("Before Total area: %f", surface_area(faces, lengths));
		
		int cnt = flipToDelaunay(faces, G, lengths);
		log("The number of flipped edges: %i", cnt);
		
		log("After Total area: %f", surface_area(faces, lengths));


		if (heatMethodAfter) {
			vector<double> distances(m.cm.vert.size());
			heatMethod(m.cm, faces, lengths, sourceVertex, distances);
			file << faces.size() << "\n";
			for (int i = 0; i < faces.size(); i++) {
				file << faces[i][0] << " " << faces[i][1] << " " << faces[i][2] << "\n";
			}
			file << m.cm.vert.size() << "\n" << sourceVertex << "\n";
			for (int i = 0; i < m.cm.vert.size(); i++)
				file << distances[i] << "\n";
			for (int i = 0; i < m.cm.vert.size(); i++)
				file << m.cm.vert[i].P().X() << " " << m.cm.vert[i].P().Y() << " " << m.cm.vert[i].P().Z() << "\n";
			file.close();
			MeshModel& nm2 = *md.addNewMesh("", "Heat Method After IDT");
			tri::Append<CMeshO, CMeshO>::MeshCopy(nm2.cm, m.cm);
			nm2.cm.face.EnableColor();

			vector<array<double, 3>> colors(distances.size());

			double minDist = *min_element(distances.begin(), distances.end());
			double maxDist = *max_element(distances.begin(), distances.end());
			double range   = maxDist - minDist;

			double radiusStep = 0.05;

			for (size_t i = 0; i < nm2.cm.face.size(); ++i) {
				CMeshO::FaceType& face = nm2.cm.face[i];

				double avgDistance = 0.0;

				for (int j = 0; j < 3; ++j) {
					int vertexIndex = face.V(j)->Index();
					avgDistance += distances[vertexIndex];
				}
				avgDistance /= 3.0;

				double normalized = (avgDistance - minDist) / range;
				int    band       = static_cast<int>(floor(normalized / radiusStep));

				double green = (band % 2 == 0) ? 255.0 : 50.0;

				face.C() = vcg::Color4b(0, static_cast<unsigned char>(green), 0, 255);
			}

			nm2.updateDataMask(MeshModel::MM_FACECOLOR);
		}

		if (showIDT) {
			MeshModel& nm3 = *md.addNewMesh("", "Heat Method Before IDT");
			for (int i = 0; i < faces.size(); i++) {
				tri::Allocator<CMeshO>::AddFace(
					nm3.cm,
					m.cm.vert[faces[i][0]].P(),
					m.cm.vert[faces[i][1]].P(),
					m.cm.vert[faces[i][2]].P());
			}
			tri::Clean<CMeshO>::RemoveDuplicateVertex(nm3.cm);
			nm3.updateBoxAndNormals();
		}

	} break;
	default: wrongActionCalled(filter);
	}
	return std::map<std::string, QVariant>();
}
MESHLAB_PLUGIN_NAME_EXPORTER(QhullPlugin)
