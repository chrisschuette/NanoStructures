/* 
 * File:   Configuration.cpp
 * Author: chris
 * 
 * Created on January 2, 2012, 2:32 PM
 */

#include "../err/exception.h"
#include "configuration.h"

#include <iostream>
#include <sstream>

using namespace config;

std::auto_ptr<Configuration> Configuration::m_instancePtr(0);

Configuration& Configuration::getInstance() {
    if(m_instancePtr.get() == 0)
    {
        m_instancePtr = std::auto_ptr<Configuration>(new Configuration);
    }
    return *m_instancePtr.get();
}

Configuration::Configuration()
{}

bool Configuration::exists(std::string path) {
    return m_configuration.exists(path.c_str());
}


void Configuration::readFile(std::string configFile) {
    m_configuration.readFile(configFile.c_str());
}

void Configuration::addSetting(std::string path, std::string name, std::string value) {
	libconfig::Setting& group = m_configuration.lookup(path);
	libconfig::Setting& setting = group.add(name, libconfig::Setting::TypeString);
	setting = value;
}

void Configuration::editSetting(std::string path, std::string name, double value) {
  std::string fullpath = "";
  if(path != "")
    fullpath = path + ".";
  fullpath += name;
  
  libconfig::Setting& setting = m_configuration.lookup(fullpath);
  setting = value;
}

double Configuration::getDouble(std::string path) {
	return m_configuration.lookup(path);
}

bool Configuration::getBool(std::string path) {
	return m_configuration.lookup(path);
}

int Configuration::getInteger(std::string path) {
	return m_configuration.lookup(path);
}

std::string Configuration::getString(std::string path){
	return m_configuration.lookup(path);
}

Configuration::~Configuration() {
}
