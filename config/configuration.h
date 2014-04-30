/* 
 * File:   Configuration.h
 * Author: chris
 *
 * Created on January 2, 2012, 2:32 PM
 */

#ifndef CONFIGURATION_H
#define	CONFIGURATION_H

#include <libconfig.h++>
#include <string>
#include <list>
#include <memory>
#include <vector>
#include <utility>

namespace config {
/*! \brief
 *
 */
class Configuration {
public:
    //Singleton
    static Configuration& getInstance();
    friend class std::auto_ptr<Configuration>;
    void readFile(std::string configFile);
    
    bool exists(std::string path);

	double getDouble(std::string path);
	bool getBool(std::string path);
	int getInteger(std::string path);
	std::string getString(std::string path);
    
    void addSetting(std::string path, std::string name, std::string value);
    void editSetting(std::string path, std::string name, double value);

    inline libconfig::Setting& lookup(const std::string path) const
   { return m_configuration.lookup(path);  }

    std::string getName() const {
        return "Configuration";
    }

private:
    //private constructor and destructor
    Configuration();
    virtual ~Configuration();
    
    static std::auto_ptr<Configuration> m_instancePtr;
     
    //config object
    libconfig::Config m_configuration;

};
}
#endif	/* CONFIGURATION_H */

