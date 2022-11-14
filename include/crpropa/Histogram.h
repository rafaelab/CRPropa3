#ifndef CRPROPA_HISTOGRAM_H
#define CRPROPA_HISTOGRAM_H


#include <string>
#include <vector>

#include <crpropa/Random.h>
#include <crpropa/Referenced.h>




namespace crpropa {

/**
 @class Histogram1D
 @brief Builds a one-dimensional histogram
 */
class Histogram1D : public Referenced {
	protected:
		std::vector<double> edges;
		std::vector<double> centres;
		std::vector<double> widths;
		std::vector<double> contents;
		std::string scale;
	public:
		Histogram1D(double vmin, double vmax, int nBins, std::string scale = "lin") {
			setScale(scale);
			initBins(vmin, vmax, nBins);
		}

		~Histogram1D(){
		}

		void initBins(double vmin, double vmax, int nBins) {
			// create empty histogram
			if (scale == "log10") {
				vmin = log10(vmin);
				vmax = log10(vmax);
				for (size_t i = 0; i < nBins + 1; i++) {
					double v = vmin + i * (vmax - vmin) / (nBins);
					edges.push_back(pow(10, v));
				}
				for (size_t i = 0; i < nBins; i++) {
					centres.push_back(pow(10, (log10(edges[i + 1]) + log10(edges[i])) / 2.));
				}
			} else {
				for (size_t i = 0; i < nBins + 1; i++) {
					double v = vmin + i * (vmax - vmin) / (nBins);
					edges.push_back(v);
				}
				for (size_t i = 0; i < nBins; i++) {
					centres.push_back((edges[i + 1] + edges[i]) / 2.);
				}
			}

			// bins start with no content
			for (size_t i = 0; i < nBins; i++) {
				widths.push_back(edges[i + 1] - edges[i]);
				contents.push_back(0);
			}
		}

		void setScale(std::string s) {
			if (s == "log")
				scale = "log10";

			if (s != "log" && s != "log10" && s != "lin")
				throw std::runtime_error("Unknown scale " + s + ".");

			scale = s;
		}

		std::string getScale() const {
			return scale;
		}

		std::vector<double> getBinEdges() const {
			return edges;
		}

		std::vector<double> getBinWidths() const {
			return widths;
		}

		std::vector<double> getBinCentres() const {
			return centres;
		}

		std::vector<double> getBinContents() const {
			return contents;
		}

		int getNumberOfBins() const {
			return contents.size();
		}

		double getSample() const {
			Random &random = Random::instance();
			size_t bin = random.randBin(contents);
			if (scale == "log") {
				return pow(10, log10(edges[bin]) + random.rand() * log10(widths[bin]));
			} else {
				return edges[bin] + random.rand() * widths[bin];
			}
		}

		void push(double v, double w = 1) {
			std::vector<double>::const_iterator it = std::lower_bound(edges.begin(), edges.end(), v);
			if (it == edges.begin() || it == edges.end())
				return;

			size_t idx = it - edges.begin();
			contents[idx] += w; 
		}

		void transformToCDF() {
			for (size_t i = 1; i < getNumberOfBins(); i++) {
				contents[i] += contents[i - 1];
			}
		}

		void clear() {
			for (size_t i = 0; i < getNumberOfBins(); i++) {
				contents[i] = 0;
			}
		}

};






// /**
//  @class Histogram1D
//  @brief Builds a one-dimensional histogram
//  */
// template <typename B, typename T>
// class Histogram1D {
// 	protected:
// 		std::vector<B> edges;
// 		std::vector<B> centres;
// 		std::vector<B> widths;
// 		std::vector<T> contents;
// 		std::string scale;
// 		T normalisation;
// 	public:
// 		Histogram1D();

// 		Histogram1D(B vmin, B vmax, int nBins, std::string scale = "lin", T normalisation = 1) {
// 			setScale(scale);
// 			setNormalisation(normalisation);
// 			initBins(vmin, vmax, nBins);
// 		}

// 		~Histogram1D(){
// 		}

// 		void initBins(B vmin, B vmax, int nBins) {
// 			// create empty histogram
// 			if (scale == "log10") {
// 				vmin = log10(vmin);
// 				vmax = log10(vmax);
// 				for (size_t i = 0; i < nBins + 1; i++) {
// 					double v = vmin + i * (vmax - vmin) / (nBins);
// 					edges.push_back(pow(10, v));
// 				}
// 				for (size_t i = 0; i < nBins; i++) {
// 					centres.push_back(pow(10, (log10(edges[i + 1]) + log10(edges[i])) / 2.));
// 				}
// 			} else {
// 				for (size_t i = 0; i < nBins + 1; i++) {
// 					double v = vmin + i * (vmax - vmin) / (nBins);
// 					edges.push_back(v);
// 				}
// 				for (size_t i = 0; i < nBins; i++) {
// 					centres.push_back((edges[i + 1] + edges[i]) / 2.);
// 				}
// 			}

// 			// bins start with no content
// 			for (size_t i = 0; i < nBins; i++) {
// 				widths.push_back(edges[i + 1] - edges[i]);
// 				contents.push_back(0);
// 			}
// 		}

// 		void setScale(std::string s) {
// 			if (s == "log")
// 				scale = "log10";

// 			if (s != "log" && s != "log10" && s != "lin")
// 				throw std::runtime_error("Unknown scale " + s + ".");

// 			scale = s;
// 		}

// 		void setNormalisation(T norm) {
// 			normalisation = norm;
// 		}

// 		T getNormalisation() const {
// 			return normalisation;
// 		}

// 		std::string getScale() const {
// 			return scale;
// 		}

// 		std::vector<B> getBinEdges() const {
// 			return edges;
// 		}

// 		std::vector<B> getBinWidths() const {
// 			return widths;
// 		}

// 		std::vector<B> getBinCentres() const {
// 			return centres;
// 		}

// 		std::vector<T> getBinContents() const {
// 			return contents;
// 		}

// 		int getNumberOfBins() const {
// 			return contents.size();
// 		}

// 		T getSample() const {
// 			Random &random = Random::instance();
// 			size_t bin = random.randBin(contents);
// 			// if (scale == "log") {
// 			// 	return pow(10, log10(edges[bin]) + random.rand() * log10(widths[bin]));
// 			// } else {
// 				return edges[bin] + random.rand() * widths[bin];
// 			// }
// 		}

// 		void push(T v, T w = 1) {
// 			std::vector<double>::const_iterator it = std::lower_bound(edges.begin(), edges.end(), v);
// 			if (it == edges.begin() || it == edges.end())
// 				return;

// 			size_t idx = it - edges.begin();
// 			contents[idx] += w; 
// 		}

// 		void transformToCDF() {
// 			for (size_t i = 1; i < getNumberOfBins(); i++) {
// 				contents[i] += contents[i - 1];
// 			}
// 		}

// 		void clear() {
// 			for (size_t i = 1; i < getNumberOfBins(); i++) {
// 				contents[i] = 0;
// 			}
// 		}

// };


} // namespace crpropa

#endif // CRPROPA_HISTOGRAM_H