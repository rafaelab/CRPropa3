#include "crpropa/Common.h"



namespace crpropa {

constexpr size_t index(size_t i, size_t j, size_t y_size) {
	return j + i * y_size;
}

std::string getDataPath(std::string filename) {
	static std::string dataPath;
	if (! dataPath.empty())
		return (std::filesystem::path(dataPath) / filename).string();

	if (const char* env_path = getenv("CRPROPA_DATA_PATH")) {
		if (std::filesystem::is_directory(env_path)) {
			dataPath = env_path;
			KISS_LOG_INFO << "getDataPath: use environment variable, " << dataPath << std::endl;
			return (std::filesystem::path(dataPath) / filename).string();
		}
	}

	#ifdef CRPROPA_INSTALL_PREFIX
		if (std::filesystem::path _path = CRPROPA_INSTALL_PREFIX "/share/crpropa"; std::filesystem::is_directory(_path)) {
			dataPath = _path.string();
			KISS_LOG_INFO << "getDataPath: use install prefix, " << dataPath << std::endl;
			return (_path / filename).string();
		}
	#endif

	if (std::filesystem::path _path = std::filesystem::current_path() / "../data"; std::filesystem::is_directory(_path)) {
		dataPath = _path.string();
		KISS_LOG_INFO << "getDataPath: use executable path, " << dataPath << std::endl;
		return (_path / filename).string();
	}

	dataPath = "data";
	KISS_LOG_INFO << "getDataPath: use default, " << dataPath << std::endl;
	return (std::filesystem::path(dataPath) / filename).string();
}

std::string getInstallPrefix() {
	std::string _path;

	#ifdef CRPROPA_INSTALL_PREFIX
		_path += CRPROPA_INSTALL_PREFIX;
	#endif

	return _path;
}

double interpolate(double x, std::span<const double> X, const std::span<const double> Y) {
	auto it = std::ranges::upper_bound(X, x);
	if (it == X.begin())
		return Y.front();
	if (it == X.end())
		return Y.back();

	size_t i = std::ranges::distance(X.begin(), it) - 1;
	return Y[i] + (x - X[i]) * (Y[i + 1] - Y[i]) / (X[i + 1] - X[i]);
}

double interpolate2d(double x, double y, std::span<const double> X, std::span<const double> Y, std::span<const double> Z) {
	if (x < X.front() or x > X.back() or y < Y.front() or y > Y.back())
		return 0;

	auto itx = std::ranges::upper_bound(X, x);
	auto ity = std::ranges::upper_bound(Y, y);

	size_t i = std::ranges::distance(X.begin(), itx) - 1;
	size_t j = std::ranges::distance(Y.begin(), ity) - 1;

	double Q11 = Z[index(i, j, Y.size())];
	double Q12 = Z[index(i, j + 1, Y.size())];
	double Q21 = Z[index(i + 1, j, Y.size())];
	double Q22 = Z[index(i + 1, j + 1, Y.size())];

	double R1 = ((X[i + 1] - x) / (X[i + 1] - X[i])) * Q11 + ((x - X[i]) / (X[i + 1] - X[i])) * Q21;
	double R2 = ((X[i + 1] - x) / (X[i + 1] - X[i])) * Q12 + ((x - X[i]) / (X[i + 1] - X[i])) * Q22;

	return ((Y[j + 1] - y) / (Y[j + 1] - Y[j])) * R1 + ((y - Y[j]) / (Y[j + 1] - Y[j])) * R2;
}

double interpolateEquidistant(double x, double lo, double hi, const  std::span<const double> Y) {
	if (x <= lo)
		return Y.front();
	if (x >= hi)
		return Y.back();

	double dx = (hi - lo) / (Y.size() - 1);
	double p = (x - lo) / dx;
	size_t i = static_cast<size_t>(p);
	
	return Y[i] + (p - i) * (Y[i + 1] - Y[i]);
}

size_t closestIndex(double x, const std::span<const double> X) {
	auto it = std::ranges::min_element(X, [x](double a, double b) {
		return std::fabs(a - x) < std::fabs(b - x);
	});
	return std::ranges::distance(X.begin(), it);
}

} // namespace crpropa

