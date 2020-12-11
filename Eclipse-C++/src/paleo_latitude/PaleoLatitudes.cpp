/*
 * PaleoLatitudes.cpp
 *
 *  Created on: 6 Jan 2017
 *      Author: Menno R.T. Fraters
 */

#include "PaleoLatitudes.h"
#include "PaleoLatitude.h"

#include <algorithm>
#include <vector>
#include <set>
#include <boost/algorithm/string.hpp>
using namespace std;
using namespace paleo_latitude;


PaleoLatitudes::PaleoLatitudes(const string& file, const bool filename) {
	_readFromFile(file, filename);
}

PaleoLatitudes::~PaleoLatitudes(){
	delete _csvdata;
	_csvdata = NULL;
}

bool PaleoLatitudes::compute(PLParameters* pl_params){

	for (const EPEntry& entry : _csvdata->getEntries()){
		// create the pl_params for this entry
		pl_params->site_latitude = entry.latitude;
		pl_params->site_longitude = entry.longitude;
		pl_params->age = entry.age;
		pl_params->age_min = entry.age_min;
		pl_params->age_max = entry.age_max;

		_pls.push_back(new PaleoLatitude(pl_params));
                bool success = false;
		try {
			success = _pls.back()->compute();
			//if (!_pls.back()->compute()){
				// Paleolatitude computation failed. Error message will have been printed to
				// the user.
				//exit(1);
			//} 
		} catch (exception& ex){
			// Something went very wrong. Print message.
			cerr << "Unexpected error computing paleolatitude: " << ex.what() << endl;
			//exit(1);
			success = false;
		}
                
		if(success == true)
		{
		  // Output is intended for human user - print the basic paleolatitude result
		  PaleoLatitude::PaleoLatitudeEntry res = _pls.back()->getPaleoLatitude();
		  const ResultEntry resultEntry(entry.sample_id,entry.latitude, entry.longitude,
				  entry.age, entry.age_min, entry.age_max,
				  res.palat,res.palat_min,res.palat_max,
		  		  res.computed_using_plate_id);
		  _entries.push_back(resultEntry);
		}
		else
		{
			double nan = std::nan("error");
            const ResultEntry resultEntry(entry.sample_id,entry.latitude, entry.longitude,
				  entry.age, entry.age_min, entry.age_max,
				  nan,nan,nan,
		  		  1001);
		  	_entries.push_back(resultEntry);
		}
	}
	return true;
}

PaleoLatitudes* PaleoLatitudes::readFromFile(const string& filename) {
	PaleoLatitudes* res = NULL;
	return res;
}

void PaleoLatitudes::_readFromFile(const string& file, const bool filename) {
	if (_csvdata != NULL) delete _csvdata;
	_csvdata = new CSVFileData<EPEntry>();
	if(filename)
		_csvdata->parseFile(file);
	else
	{
		string text = file;
		boost::replace_all(text,"\\n","\n");
		istringstream _istringstream (text);
		_csvdata->parseStream(&_istringstream);
	}
}


vector<unsigned int> PaleoLatitudes::getRelevantAges(const PLPlate* plate, unsigned int min, unsigned int max){
	vector<unsigned int> res;

	int left_outside_age = -1;
	int right_outside_age = -1;

	for (const EPEntry& entry : _csvdata->getEntries()){

		if (entry.age >= min && entry.age <= max){
			if (!res.empty() && res.back() == entry.age){
				// Duplicate age entry in CSV (possible for cross-over between relative plates)
				continue;
			}

			// CSV entry falls in between the min/max: this age should be included
			// in the ages.
			res.push_back(entry.age);
		}

		if (entry.age <= min){
			// Entry is on window border, or just outside window. Let's
			// see whether we need it.
			if (left_outside_age == -1 || (unsigned int) left_outside_age < entry.age){
				// The previously found left-outside entry is less relevant than
				// this entry
				left_outside_age = entry.age;
			}
		}

		if (entry.age >= max){
			// Entry is on/outside window border on right
			if (right_outside_age == -1 || (unsigned int) right_outside_age > entry.age){
				right_outside_age = entry.age;
			}
		}
	}

	if (left_outside_age >= 0 && (unsigned int) left_outside_age < min){
		// Include left-outside
		res.push_back(left_outside_age);
	} // else: no left-outside, or left outside is on window border (and therefore already included)


	if (right_outside_age >= 0 && (unsigned int) right_outside_age > max){
		// Include left-outside
		res.push_back(right_outside_age);
	} // else: no right-outside, or is on window border (and therefore already included)

	// Make sure the years are sorted from small (recent) to large (less recent)
	sort(res.begin(), res.end());

	return res;
}


void PaleoLatitudes::EPEntry::set(unsigned int col_index, const string value, unsigned int lineno){
	// Column index 0: sample ID (string)
	// Column index 1: latitude (double)
	// Column index 2: longitude (double
	// Column index 3: age (uint, in Myr)
	// Column index 4: age min (uint, in Myr)
	// Column index 5: age max (uint, in Myr)

	if (col_index == 0) sample_id = value;
	if (col_index == 1) this->parseString(value, this->latitude);
	if (col_index == 2) this->parseString(value, this->longitude);
	if (col_index == 3) this->parseString(value, this->age);
	if (col_index == 4) this->parseString(value, this->age_min);
	if (col_index == 5) this->parseString(value, this->age_max);
}

size_t PaleoLatitudes::EPEntry::numColumns() const {
	return 6;
}


bool PaleoLatitudes::EPEntry::compareByAge(const EPEntry* a, const EPEntry* b) {
	return (a->age < b->age);
}


vector<const PaleoLatitudes::EPEntry*> PaleoLatitudes::getEntries(unsigned int plate_id) const {
	vector<const PaleoLatitudes::EPEntry*> res;

	for (const EPEntry& entry : _csvdata->getEntries()){
		res.push_back(&entry);
	}
	return res;
}

vector<PaleoLatitudes::EPEntry*> PaleoLatitudes::getEntries(unsigned int plate_id, unsigned int age) const {
	vector<PaleoLatitudes::EPEntry*> res;

	if (res.size() == 0){
		Exception ex;
		ex << "No entry for age=" << age << " and plate_id=" << plate_id << " found in Euler pole table?";
		throw ex;
	}

	return res;
}


const vector<PaleoLatitudes::EPEntry>& PaleoLatitudes::getAllEntries() const {
	return _csvdata->getEntries();
}

void PaleoLatitudes::writeCSV(const string& filename) {
	ofstream out(filename);
	writeCSV(out);
	out.close();
}

void PaleoLatitudes::writeCSV(ostream& output_stream) {
	_requireResult();
	output_stream << "# sample id, latitude, longitude, age, age min, age_max, palat, palat min, palat max, plate id" << endl;
	output_stream.setf(ios::fixed, ios::floatfield);

	for (const ResultEntry entry : _entries){
		output_stream << entry.sample_id << ",";

		output_stream.precision(5);
		output_stream << entry.latitude << ",";

		output_stream.precision(5);
		output_stream << entry.longitude << ",";

		output_stream.precision(5);
		output_stream << entry.age << ",";

		output_stream.precision(5);
		output_stream << entry.age_min << ",";

		output_stream.precision(5);
		output_stream << entry.age_max << ",";

		output_stream.precision(5);
		output_stream << entry.palat << ",";

		output_stream.precision(5);
		output_stream << entry.palat_min << ",";

		output_stream.precision(5);
		output_stream << entry.palat_max << ",";

		output_stream.precision(5);
		output_stream << entry.computed_using_plate_id;

		output_stream << endl;
	}
}

void PaleoLatitudes::_requireResult() const {
	if (_entries.size() == 0) throw Exception("PaleoLatitude not computed - call compute() first");
}

