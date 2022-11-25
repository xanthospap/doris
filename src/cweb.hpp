#ifndef __DSO_CPP_WEB_INTERFACE__
#define __DSO_CPP_WEB_INTERFACE__

namespace dso {
/// @brief Download remote file
/// @param[in] url The url, remote file to download
/// @param[in] local The local filename to save remote file to
/// @param[in] force Force download if local already exists
/// @return An integer denoting:
///        int < 0 : file already exists and force was false; aka we did not
///                  download the file but the file already exists
///              0 : file downloaded
///            > 0 : error
int http_get(const char *url, const char *local, bool force = false) noexcept;
} // namespace dso
#endif
