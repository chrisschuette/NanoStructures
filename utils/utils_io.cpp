#include "utils_io.h"
#include <iostream>
#include <sys/stat.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <cstdarg>


namespace utils {
namespace io {
void createDirectoryStructure(std::string filePath) {
    unsigned found = filePath.find_last_of('/');
    std::string filename;
    std::string directory;
    if(found >= filePath.length()) {
        directory = "";
        filename = filePath;
    } else {
        directory = filePath.substr(0,found);

        filename = filePath.substr(found+1);
    }
    if(directory.size() > 0) {
        std::cout << "Creating " << directory << std::endl;
        mkpath(directory.c_str(), 0755);
    }

}
typedef struct stat Stat;

static int do_mkdir(const char *path, mode_t mode)
{
    Stat            st;
    int             status = 0;

    if (stat(path, &st) != 0)
    {
        /* Directory does not exist. EEXIST for race condition */
        if (mkdir(path, mode) != 0 && errno != EEXIST)
            status = -1;
    }
    else if (!S_ISDIR(st.st_mode))
    {
        errno = ENOTDIR;
        status = -1;
    }

    return(status);
}

/**
** mkpath - ensure all directories in path exist
** Algorithm takes the pessimistic view and works top-down to ensure
** each directory in path exists, rather than optimistically creating
** the last element and working backwards.
*/
int mkpath(const char *path, mode_t mode)
{
    char           *pp;
    char           *sp;
    int             status;
    char           *copypath = strdup(path);

    status = 0;
    pp = copypath;
    while (status == 0 && (sp = strchr(pp, '/')) != 0)
    {
        if (sp != pp)
        {
            /* Neither root nor double slash in path */
            *sp = '\0';
            status = do_mkdir(copypath, mode);
            *sp = '/';
        }
        pp = sp + 1;
    }
    if (status == 0)
        status = do_mkdir(path, mode);
    free(copypath);
    return (status);
}

void appendToFile(std::string filename, int num, ...) {
    std::ofstream file(filename.c_str(), std::ios::app);
    va_list arguments;                     // A place to store the list of arguments
    va_start ( arguments, num );           // Initializing arguments to store all values after num
    for ( int x = 0; x < num-1; x++ ) {
        file << va_arg ( arguments, double ) << " "; // Adds the next value in argument list to sum.
    }
    file << va_arg ( arguments, double ) << std::endl;
    va_end ( arguments );                  // Cleans up the list
    file.close();
}


}
}
