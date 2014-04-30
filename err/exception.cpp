#include "exception.h"

namespace error {
Exception::Exception(std::string message) : std::exception(), m_message(message)
{
}
}
