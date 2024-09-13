#ifndef ReadCorrection_h
#define ReadCorrection_h

#include <string>
#include <vector>
#include <stdexcept>
#include <cctype>
#include <cstdint>  // for uint32_t
#include <map>

struct Correction
{
	size_t startIndex;
	size_t endIndex;
	std::string corrected;
	std::string cigar;
};

std::string getCorrected(const std::string& raw, const std::vector<Correction>& corrections, size_t maxOverlap);

std::map<size_t, size_t> get_aligned_positions(
    size_t query_start, size_t reference_start, const std::string& cigar
);

#endif
