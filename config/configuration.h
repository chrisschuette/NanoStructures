/* 
 * File:   Configuration.h
 * Author: chris
 *
 * Created on January 2, 2012, 2:32 PM
 */

#ifndef CONFIGURATION_H
#define	CONFIGURATION_H

/**
 * @file
 *
 * @ingroup config
 *
 * @brief A configuration object which is populated from a configuration file
 *        and allows to query and modify configuration values.
 * 
 * Essentially this class is just a front-end for the libconfig++ library. For 
 * convenience the configuration is a singleton object.
 */

/*
 * Copyright (c) 2014 Christoph Schuette.
 *
 * The license and distribution terms for this file may be
 * found in the file LICENSE in this distribution
 */

#include <libconfig.h++>
#include <string>
#include <list>
#include <memory>
#include <vector>
#include <utility>

namespace config {

    /**
     * @brief A configuration object which is populated from a configuration file
     *        and allows to query and modify configuration values.
     * 
     * Essentially this class is just a front-end for the libconfig++ library. For 
     * convenience the configuration is a singleton object.
     */
    class Configuration {
    public:
        /**
         * \brief returns a reference to the only configuration instance
         *
         * If the configuration has not been accessed previously an instance is created.
         * The auto_ptr ensures that the object is properly destructed when the 
         * the application closes.
         * 
         * @return reference to the only configuration instance 
         */
        static Configuration& getInstance();
        friend class std::auto_ptr<Configuration>;

        /**
         * \brief parses configuration file configFile and loads the data into
         *        object.
         * The configuration file needs to written in the libconfig++ format,
         * see <a href="">libconfig++</a>.
         * 
         * @param[in] configFile path to the configuration file
         */
        void readFile(std::string configFile);

        /**
         * \brief determines whether a given configuration value exists.
         * 
         * Libconfig++ configuration files have the structure of a tree, as they
         * support groups of configuration values which themselves can be elements.
         * A path to a configuration value within this tree uses a . as the separator.
         * 
         * example: "group1.nestedgroupA.value"
         * 
         * @param[in] path path to the configuration value
         * @return true if the value is defined otherwise false.
         */
        bool exists(std::string path);

        /**
         * \brief retrieves the configuration value found at a given configuration path.
         * 
         * If the lookup fails the function throws a SettingNotFoundException
         * exception.
         * 
         * @param[in] path path to the configuration value
         * @return double configuration value
         */
        double getDouble(std::string path);
        
        /**
         * \brief retrieves the configuration value found at a given configuration path.
         * 
         * If the lookup fails the function throws a SettingNotFoundException
         * exception.
         * 
         * @param[in] path path to the configuration value
         * @return bool configuration value
         */
        bool getBool(std::string path);
        
        /**
         * \brief retrieves the configuration value found at a given configuration path.
         * 
         * If the lookup fails the function throws a SettingNotFoundException
         * exception.
         * 
         * @param[in] path path to the configuration value
         * @return integer configuration value
         */
        int getInteger(std::string path);
        
        /**
         * \brief retrieves the configuration value found at a given configuration path.
         * 
         * If the lookup fails the function throws a SettingNotFoundException
         * exception.
         * 
         * @param[in] path path to the configuration value
         * @return string configuration value
         */
        std::string getString(std::string path);

        /**
         * \brief adds the string setting named <i>name</i> to the configuration tree at
         *        <i>path</i> with the value <i>value</i>
         * 
         * @param[in] path path to the configuration value
         * @param[in] name  name of the setting
         * @param[in] value value of the setting
         */
        void addSetting(std::string path, std::string name, std::string value);
        
        /**
         * \brief modifies the string setting named <i>name</i> to the configuration tree at
         *        <i>path</i> with the value <i>value</i>
         * 
         * @param[in] path path to the configuration value
         * @param[in] name  name of the setting
         * @param[in] value value of the setting
         */
        void editSetting(std::string path, std::string name, double value);

        /**
         * \brief retrieves the configuration value found at a given configuration path
         *        as a reference to a libconfig::Setting object.
         * libconfig::Setting objects support automatic conversion to most common
         * data types.
         * 
         * @param[in] path path to the configuration value
         * @return reference to a Setting object
         */
        inline libconfig::Setting& lookup(const std::string path) const {
            return m_configuration.lookup(path);
        }

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

