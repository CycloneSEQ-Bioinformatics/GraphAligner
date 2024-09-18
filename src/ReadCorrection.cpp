#include "ThreadReadAssertion.h"
#include "ReadCorrection.h"

std::string toUpper(std::string seq)
{
	for (auto& c : seq)
	{
		c = toupper(c);
	}
	return seq;
}

std::string toLower(std::string seq)
{
	for (auto& c : seq)
	{
		c = tolower(c);
	}
	return seq;
}

size_t getLongestOverlap(const std::string& left, const std::string& right, size_t maxOverlap)
{
	if (left.size() < maxOverlap) maxOverlap = left.size();
	if (right.size() < maxOverlap) maxOverlap = right.size();
	for (size_t i = maxOverlap; i > 0; i--)
	{
		bool match = true;
		for (size_t a = 0; a < i && match; a++)
		{
			if (left[left.size() - maxOverlap + a] != right[a]) match = false;
		}
		if (match) return i;
	}
	return 0;
}

std::string getCorrected(const std::string& raw, const std::vector<Correction>& corrections, size_t maxOverlap)
{
	std::string result;
	size_t currentEnd = 0;
	for (size_t i = 0; i < corrections.size(); i++)
	{
		assert(i == 0 || corrections[i].startIndex >= corrections[i-1].startIndex);
		if (corrections[i].startIndex < currentEnd)
		{
			if(corrections[i].endIndex <= currentEnd)continue;
			/*
			if(currentEnd - corrections[i].startIndex > maxOverlap){
				std::string corrected_str = corrections[i].corrected.substr(currentEnd - corrections[i].startIndex - maxOverlap);
				size_t overlap = getLongestOverlap(result, corrected_str, maxOverlap);
				if (readName == "WTB4000011D2-202307071450251_20230707193002_00764_0043826830_13.79"){
					std::cout<<"overlap:\t"<<overlap<<"\t"<<corrected_str<<std::endl;
				}
				result += toUpper(corrected_str.substr(overlap));
			}else{
				size_t overlap = getLongestOverlap(result, corrections[i].corrected, maxOverlap);
				if (readName == "WTB4000011D2-202307071450251_20230707193002_00764_0043826830_13.79"){
					std::cout<<"overlap:\t"<<overlap<<std::endl;
				}
				result += toUpper(corrections[i].corrected.substr(overlap));
			}
			*/
			std::map mapping = get_aligned_positions(corrections[i].startIndex, 0, corrections[i].cigar); // mapping: raw -> corrected
			if(mapping.rbegin()->first <= currentEnd)continue;
			size_t count_bad_base = 0;
			while(mapping.find(currentEnd-count_bad_base) == mapping.end()){
				count_bad_base ++;
			}
			result = result.substr(0, result.length()-count_bad_base);
			result += toUpper(corrections[i].corrected.substr(mapping.at(currentEnd-count_bad_base)));
		}
		else if (corrections[i].startIndex > currentEnd)
		{
			result += toLower(raw.substr(currentEnd, corrections[i].startIndex - currentEnd));
			result += toUpper(corrections[i].corrected);
		}
		else
		{
			assert(corrections[i].startIndex == currentEnd);
			result += toUpper(corrections[i].corrected);
		}
		currentEnd = corrections[i].endIndex;
	}
	if (currentEnd < raw.size()) result += toLower(raw.substr(currentEnd));
	return result;
}


std::map<size_t, size_t> get_aligned_positions(
    size_t query_start, size_t reference_start, const std::string& cigar
) {
    size_t query_position = query_start;
    size_t reference_position = reference_start;
    std::map<size_t, size_t> mapping;  // ref pos -> qry pos

    uint32_t length = 0;
    char op = '\0';

    for (size_t i = 0; i < cigar.length(); ++i) {
        if (isdigit(cigar[i])) {
            length = length * 10 + (cigar[i] - '0');
        } else {
            op = cigar[i];
            switch (op) {
                case 'M':  // Alignment match (can be a sequence match or mismatch)
                case '=':  // Sequence match
                case 'X':  // Sequence mismatch
                    for (uint32_t j = 0; j < length; ++j) {
                        mapping[query_position] = reference_position;
                        query_position++;
                        reference_position++;
                    }
                    break;

                case 'I':  // Insertion to the reference
                    query_position += length;
                    break;
                
                case 'D':  // Deletion from the reference
                    for (uint32_t j = 0; j < length; ++j) {
                        mapping[query_position] = reference_position;
                        reference_position++;
                    }
                    break;
                
                case 'N':  // Skipped region from the reference
                case 'S':  // Soft clipping
                    query_position += length;
                    break;
                
                case 'H':  // Hard clipping (clipped sequences not present in the alignment)
                    break;
                
                case 'P':  // Padding (silent deletion from padded reference)
                    reference_position += length;
                    break;
                
                default:
                    throw std::invalid_argument("Unsupported CIGAR operation: " + std::string(1, op));
            }
            // 重置长度，准备解析下一个操作
            length = 0;
        }
    }

    return mapping;
}