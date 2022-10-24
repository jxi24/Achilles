#ifndef ACHILLES_EXCEPTION_HH
#define ACHILLES_EXCEPTION_HH

#include <stdexcept>
#include "fmt/format.h"

namespace achilles {

    struct achilles_error : public std::exception {
        std::string msg;
        achilles_error() : msg{} {}
        achilles_error(const std::string &_msg) : msg{_msg} {}
        const char *what() const noexcept { return msg.c_str(); }
    };

    struct option_error : public achilles_error {
        option_error(const std::string &err_msg) {
            msg = fmt::format("Achilles::OptionError: {}", err_msg); 
        }
    };

    struct parser_error : public achilles_error {
        parser_error(const std::string &err_msg) {
            msg = fmt::format("Achilles::ParserError: {}", err_msg); 
        }
    };

    struct not_implemented_error : public achilles_error {
        not_implemented_error(const std::string &err_msg) {
            msg = fmt::format("Achilles::NotImplementedError: {}", err_msg); 
        }
    };

    struct cascade_error : public achilles_error {
        cascade_error(const std::string &err_msg) {
            msg = fmt::format("Achilles::CascadeError: {}", err_msg);
        }
    };

    struct beam_error : public achilles_error {
        beam_error(const std::string &err_msg) {
            msg = fmt::format("Achilles::BeamSpectrumError: {}", err_msg);
        }
    };

    struct nuclear_error : public achilles_error {
        nuclear_error(const std::string &model, const std::string &err_msg) {
            msg = fmt::format("Achilles::NuclearError:{}: {}", model, err_msg);
        }
    };

    struct nucleus_error : public achilles_error {
        nucleus_error(const std::string &err_msg) {
            msg = fmt::format("Achilles::NucleusError: {}", err_msg);
        }
    };

}


#endif
