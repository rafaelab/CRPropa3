#pragma once

#include <cmath>
#include <numeric>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>
#include <algorithm>

#include "crpropa/Random.h"
#include "crpropa/Referenced.h"




namespace crpropa {


/**
 * \addtogroup Tools
 * @{
 */



/**
 * @class Histogram1D
 * @brief Builds a one-dimensional histogram.
 * @details The histogram is defined by its edges, centres, widths and contents.
 * The scale can be set to "lin" for linear or "log10" for logarithmic scale. 
 * The histogram can be sampled, normalised, integrated, transformed to PDF or CDF
 */
class Histogram1D : public Referenced {
	protected:
		std::vector<double> edges;
		std::vector<double> centres;
		std::vector<double> widths;
		std::vector<double> contents;
		std::string scale;
	public:
		Histogram1D(double vmin, double vmax, int nBins, std::string_view scale = "lin")
			: scale(scale) {
			initBins(vmin, vmax, nBins);
		}

		~Histogram1D() = default;

		void initBins(double vmin, double vmax, int nBins) {
			edges.resize(nBins + 1);
			centres.resize(nBins);
			widths.resize(nBins);
			contents.resize(nBins, 0);

			if (scale == "log10") {
				vmin = log10(vmin);
				vmax = log10(vmax);
				for (int i = 0; i <= nBins; ++i) {
					edges[i] = pow(10, vmin + i * (vmax - vmin) / nBins);
				}
				for (int i = 0; i < nBins; ++i) {
					centres[i] = pow(10, (log10(edges[i + 1]) + log10(edges[i])) / 2.);
					widths[i] = edges[i + 1] - edges[i];
				}
			} else {
				for (int i = 0; i <= nBins; ++i) {
					edges[i] = vmin + i * (vmax - vmin) / nBins;
				}
				for (int i = 0; i < nBins; ++i) {
					centres[i] = (edges[i + 1] + edges[i]) / 2.;
					widths[i] = edges[i + 1] - edges[i];
				}
			}
		}

		constexpr void setScale(std::string_view s) {
			if (s == "log") {
				scale = "log10";
			} else if (s == "log10" || s == "lin") {
				scale = s;
			} else {
				throw std::runtime_error("Unknown scale " + std::string(s) + ".");
			}
		}

		[[nodiscard]] std::string getScale() const {
			return scale;
		}

		[[nodiscard]] std::vector<double> getBinEdges() const {
			return edges;
		}

		[[nodiscard]] std::vector<double> getBinWidths() const {
			return widths;
		}

		[[nodiscard]] std::vector<double> getBinCentres() const {
			return centres;
		}

		[[nodiscard]] std::vector<double> getBinContents() const {
			return contents;
		}

		[[nodiscard]] int getNumberOfBins() const {
			return static_cast<int>(contents.size());
		}

		[[nodiscard]] double getSample() const {
			Random& random = Random::instance();
			size_t bin = random.randBin(contents);
			if (scale == "log") {
				return pow(10, log10(edges[bin]) + random.rand() * log10(widths[bin]));
			} else {
				return edges[bin] + random.rand() * widths[bin];
			}
		}

		void push(double v, double w = 1) {
			auto it = std::ranges::lower_bound(edges, v);
			if (it == edges.begin() || it == edges.end())
				return;

			size_t idx = std::distance(edges.begin(), it) - 1;
			contents[idx] += w;
		}

		void normalise(double norm) {
			for (auto& content : contents) {
				content /= norm;
			}
		}

		void normalize(double norm) {
			normalise(norm);
		}

		[[nodiscard]] double sum() const {
			return std::accumulate(contents.begin(), contents.end(), 0.0);
		}

		[[nodiscard]] double integrate() const {
			if (scale == "log") {
				return std::transform_reduce(contents.begin(), contents.end(), widths.begin(), 0.0, std::plus<>(), [](double c, double w) {
					return c / w * log(10);
				});
			} else {
				return sum();
			}
		}

		void transformToPDF() {
			double integral = integrate();
			for (auto& content : contents) {
				content /= integral;
			}
		}

		void transformToCDF() {
			std::partial_sum(contents.begin(), contents.end(), contents.begin());
		}

		void clear() {
			std::fill(contents.begin(), contents.end(), 0);
		}
};


/** @} 
 */

} // namespace crpropa
