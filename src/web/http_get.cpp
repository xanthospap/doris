#include "cweb.hpp"
#include <curl/curl.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

/// see https://curl.se/libcurl/c/url2file.html

static size_t write_data(void *ptr, size_t size, size_t nmemb, void *stream) {
  size_t written = fwrite(ptr, size, nmemb, (FILE *)stream);
  return written;
}

int dso::http_get(const char *url, const char *local) noexcept {
  CURL *curl_handle;
  FILE *pagefile;
  int status = 1;

  curl_global_init(CURL_GLOBAL_ALL);

  //  init the curl session
  curl_handle = curl_easy_init();

  //  set URL to get here
  curl_easy_setopt(curl_handle, CURLOPT_URL, url);

  // Switch on full protocol/debug output while testing
  // curl_easy_setopt(curl_handle, CURLOPT_VERBOSE, 1L);

  //  request failure on HTTP response >= 400
  curl_easy_setopt(curl_handle, CURLOPT_FAILONERROR, 1L);

  // disable progress meter, set to 0L to enable it
  curl_easy_setopt(curl_handle, CURLOPT_NOPROGRESS, 1L);

  //  send all data to this function
  curl_easy_setopt(curl_handle, CURLOPT_WRITEFUNCTION, write_data);

  //  open the file
  pagefile = fopen(local, "wb");
  if (pagefile) {

    //  write the page body to this file handle
    curl_easy_setopt(curl_handle, CURLOPT_WRITEDATA, pagefile);

    //  get it!
    status = curl_easy_perform(curl_handle);

    //  close the header file
    fclose(pagefile);
  }

  //  cleanup curl stuff
  curl_easy_cleanup(curl_handle);

  curl_global_cleanup();

  return status;
}
