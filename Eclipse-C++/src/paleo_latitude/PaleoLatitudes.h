/*
 * PaleoLatitudes.h
 *
 *  Created on: 6 Jan 2017
 *      Author: Menno R.T. Fraters
 */

#ifndef PALEOLATITUDES_H_
#define PALEOLATITUDES_H_

#include "PaleoLatitude.h"
#include "PLParameters.h"
#include <set>
#include "../util/CSVFileData.h"

namespace paleo_latitude {

class PLPlate;

/**
 * Data structure that processes and holds Euler pole reconstructions from a specific
 * model. It reads data from a single CSV file (e.g. 'euler-torsvik-2012.csv'), and
 * expects the following columns in this order:
 *  - sample ID (std::sting)
 *  - latitude (double)
 *  - longitude (double)
 *  - age (uint, in Myr)
 *  - age min (uint, in Myr)
 *  - age max (uint, in Myr)
 *
 */
class PaleoLatitudes {
public:
	struct ResultEntry
	{
		ResultEntry(std::string sample_id, double latitude, double longitude,
		unsigned int age, unsigned int age_min, unsigned int age_max,
		double palat, double palat_min, double palat_max, 
		unsigned int computed_using_plate_id) 
		: sample_id(sample_id), latitude(latitude), longitude(longitude), 
				age(age), age_min(age_min), age_max(age_max), 
				palat(palat), palat_min(palat_min), palat_max(palat_max), 
				computed_using_plate_id(computed_using_plate_id){};
		std::string sample_id = "";
		double latitude = 0;
		double longitude = 0;
		unsigned int age = 0;
		unsigned int age_min = 0;
		unsigned int age_max = 0;
		double palat = 0;
		double palat_min = 0;
		double palat_max = 0;
		unsigned int computed_using_plate_id = 0;
	};
	
	struct EPEntry : public CSVFileData<EPEntry>::StringEntry {
		EPEntry(CSVFileData<EPEntry>* parent, unsigned int line_no) : CSVFileData<EPEntry>::StringEntry(parent, line_no){}
		void set(unsigned int col_index, const string value, unsigned int lineno) override;
		size_t numColumns() const override;

		static bool compareByAge(const EPEntry* a, const EPEntry* b);

		std::string sample_id = "";
		double latitude = 0;
		double longitude = 0;
		unsigned int age = 0;
		unsigned int age_min = 0;
		unsigned int age_max = 0;

	};

	PaleoLatitudes(const string& file, const bool filename);
	virtual ~PaleoLatitudes();
	
	/**
	 * Equivalent to the compute of PaleoLatitude
	 */
	bool compute(PLParameters* pl_params);

	/**
	 * Returns all ages relevant to a paleolatitude query from min to max.
	 */
	vector<unsigned int> getRelevantAges(const PLPlate* plate, unsigned int min, unsigned int max);

	/**
	 * Returns the Euler poles for a given plate and age. Often, the result will be a single
	 * EPEntry, but in rare cases two entries will be returned. This happens at the cross-over
	 * point at which rotation is expressed relative to two plates (e.g. plate 102 at 320 Ma in Torsvik).
	 */
	vector<EPEntry*> getEntries(unsigned int plate_id, unsigned int age) const;

	/**
	 * Returns all Euler poles for a given plate ID.
	 */
	vector<const EPEntry*> getEntries(unsigned int plate_id) const;

	/**
	 * Returns all plate IDs found in the Euler data
	 */
	set<unsigned int> getPlateIds() const;

	//vector<const EPEntry*> getEntries(const PLPlate* plate, unsigned int age) const;

	const vector<EPEntry>& getAllEntries() const;
	
	void writeCSV(const string& filename);
	void writeCSV(ostream& output_stream);

	static PaleoLatitudes* readFromFile(const string& filename);

private:
	void _requireResult() const;

	void _readFromFile(const string& file, const bool filename);

	CSVFileData<EPEntry>* _csvdata = NULL;
	
	vector<PaleoLatitude*> _pls;
	vector<ResultEntry> _entries;
};

};

#endif /* PALEOLATITUDES_H_ */
