/// @file
/// @brief すぐ使える実用的な簡易ロガー
/// @description LOGI, LOGD, LOGW, LOGE の4つのマクロが定義され、それぞれをログ出力ストリームと見做して LOGI << "hoge"; のように使う。
/// LOGI ... info  レベルのログメッセージ出力（白色）
/// LOGD ... debug レベルのログメッセージ出力（緑色）
/// LOGW ... warn  レベルのログメッセージ出力（黄色）
/// LOGE ... error レベルのログメッセージ出力（赤色）
/// @note
/// (a): USE_UTIL_LOG_EASY_LOGGER_BOOST_CPU_TIMER が定義されている場合、時間計測に boost::timer::cpu_timer を使用する（ boost.timer 関連のライブラリーのリンクが必要 ）
/// (b): USE_UTIL_LOG_EASY_LOGGER_BOOST_CHRONO が定義されている場合、時間計測に boost::chrono::steady_clock を使用する ( boost.chrono 関連のライブラリーのリンクが必要)
/// (c): USE_UTIL_LOG_EASY_LOGGER_STD_CHRONO が定義されている場合、時間計測に std::chrono::steady_clock を使用する（ 外部ライブラリーのリンクは不要だが、処理系によっては分解能が不足する ）
/// (d): (a), (b), (c) の何れも定義されていない場合、時間計測に util::chrono::default_clock を使用する（ 外部ライブラリーは不要、windows処理系でもQueryPerformanceCounterを内部的に使用する ）
/// (f): DISABLE_UTIL_LOG_EASY_LOGGER が定義されている場合、全てのログ出力は事実上無効になる。
/// UTIL_LOG_EASY_LOGGER_GET_{ FUNCTION | FILE | LINE } を事前に定義しておくとユーザー定義の何かそういうものに置き換え可能。（ "" をユーザー定義すれば出力から消す事も可能。 ）
/// UTIL_LOG_EASY_LOGGER_{INFO|DEBUG|WARN|ERROR}_PREFIX を事前に定義しておくと [ info ] 的な部分をユーザー定義に置き換え可能。（ "" をユーザー定義すれば出力から消す事も可能。 ）
//https://qiita.com/usagi/items/d4aec8d3f748f4ba9d6a
//https://github.com/usagi/usagi/blob/master/include/usagi/log/easy_logger.hxx

#pragma once

#ifdef DISABLE_UTIL_LOG_EASY_LOGGER


#define LOGI ::util::log::easy_logger::log_null()
#define LOGD ::util::log::easy_logger::log_null()
#define LOGW ::util::log::easy_logger::log_null()
#define LOGE ::util::log::easy_logger::log_null()

#else

#define __FILENAME__ (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)

#ifndef UTIL_LOG_EASY_LOGGER_GET_FILE
  #define UTIL_LOG_EASY_LOGGER_GET_FILE __FILENAME__
#endif

#ifndef UTIL_LOG_EASY_LOGGER_GET_LINE
  #define UTIL_LOG_EASY_LOGGER_GET_LINE __LINE__
#endif

#ifndef UTIL_LOG_EASY_LOGGER_GET_FUNCTION
  #if defined( __Clang__ ) || defined( __GNUC__ )
    #define UTIL_LOG_EASY_LOGGER_GET_FUNCTION __FUNCTION__
  #else
    #define UTIL_LOG_EASY_LOGGER_GET_FUNCTION __func__
  #endif
#endif

#ifndef UTIL_LOG_EASY_LOGGER_INFO_PREFIX
  #define UTIL_LOG_EASY_LOGGER_INFO_PREFIX  " [ INFO  ] "
#endif

#ifndef UTIL_LOG_EASY_LOGGER_DEBUG_PREFIX
  #define UTIL_LOG_EASY_LOGGER_DEBUG_PREFIX " [ DEBUG ] "
#endif

#ifndef UTIL_LOG_EASY_LOGGER_WARN_PREFIX
  #define UTIL_LOG_EASY_LOGGER_WARN_PREFIX  " [ WARN ] "
#endif

#ifndef UTIL_LOG_EASY_LOGGER_ERROR_PREFIX
  #define UTIL_LOG_EASY_LOGGER_ERROR_PREFIX " [ ERROR ] "
#endif


#define LOGI ::util::log::easy_logger::log_intermediate::make_info \
  ( UTIL_LOG_EASY_LOGGER_GET_FILE \
  , UTIL_LOG_EASY_LOGGER_GET_LINE \
  , UTIL_LOG_EASY_LOGGER_GET_FUNCTION \
  ) << UTIL_LOG_EASY_LOGGER_INFO_PREFIX
#define LOGD ::util::log::easy_logger::log_intermediate::make_debug\
  ( UTIL_LOG_EASY_LOGGER_GET_FILE \
  , UTIL_LOG_EASY_LOGGER_GET_LINE \
  , UTIL_LOG_EASY_LOGGER_GET_FUNCTION \
  ) << UTIL_LOG_EASY_LOGGER_DEBUG_PREFIX
#define LOGW ::util::log::easy_logger::log_intermediate::make_warn \
  ( UTIL_LOG_EASY_LOGGER_GET_FILE \
  , UTIL_LOG_EASY_LOGGER_GET_LINE \
  , UTIL_LOG_EASY_LOGGER_GET_FUNCTION \
  ) << UTIL_LOG_EASY_LOGGER_WARN_PREFIX
#define LOGE ::util::log::easy_logger::log_intermediate::make_error\
  ( UTIL_LOG_EASY_LOGGER_GET_FILE \
  , UTIL_LOG_EASY_LOGGER_GET_LINE \
  , UTIL_LOG_EASY_LOGGER_GET_FUNCTION \
  ) << UTIL_LOG_EASY_LOGGER_ERROR_PREFIX

#endif

#ifndef UTIL_LOG_EASY_LOGGER_OUT
  #define UTIL_LOG_EASY_LOGGER_OUT ::std::cout
#endif

#ifndef UTIL_LOG_EASY_LOGGER_FLUSH
  #define UTIL_LOG_EASY_LOGGER_FLUSH ::std::flush
#endif

#include <string>
#include <iostream>
#include <iomanip>

namespace util::log::easy_logger
{
  using namespace std;

  
#ifdef DISABLE_UTIL_LOG_EASY_LOGGER
  
  struct log_null
  {
    template < typename T >
    decltype( auto ) operator<<( const T& ) { return *this; }
  };
  
#else
  
  
  struct log_intermediate
  {
    stringstream buffer;
    
    const char* prefix;
    const char* suffix;
    
    const char* source;
    const size_t line;
    const char* function;
    
    static constexpr auto prefix_i = "";
    static constexpr auto prefix_d = "";
    static constexpr auto prefix_w = "";
    static constexpr auto prefix_e = "";
    
    static constexpr auto suffix_i = "";
    static constexpr auto suffix_d = "";
    static constexpr auto suffix_w = "";
    static constexpr auto suffix_e = "";
    
    static constexpr auto log_time_width = 16;
    
    log_intermediate( log_intermediate&& a )
      : buffer( move( a.buffer ) )
      , prefix( a.prefix )
      , suffix( a.suffix )
      , source( a.source )
      , line( a.line )
      , function( a.function )
    { }
    
    log_intermediate( const char* p, const char* s, const char* source_, const size_t line_, const char* function_ )
      : prefix( p )
      , suffix( s )
      , source( source_ )
      , line( line_ )
      , function( function_ )
    { }
    
    static auto make_info ( const char* s, const size_t l, const char* f ) { return log_intermediate( prefix_i, suffix_i, s, l, f ); }
    static auto make_debug( const char* s, const size_t l, const char* f ) { return log_intermediate( prefix_d, suffix_d, s, l, f ); }
    static auto make_warn ( const char* s, const size_t l, const char* f ) { return log_intermediate( prefix_w, suffix_w, s, l, f ); }
    static auto make_error( const char* s, const size_t l, const char* f ) { return log_intermediate( prefix_e, suffix_e, s, l, f ); }
    
    template < typename T >
    decltype( auto ) operator<<( const T& in ) noexcept
    {
      try
      { return buffer << in; }
      catch ( const exception& e )
      { return buffer << "<<<<<exception on " << __PRETTY_FUNCTION__ << " what=" << e.what() << ">>>>>"; }
      catch ( ... )
      { return buffer << "<<<<<exception on " << __PRETTY_FUNCTION__ << " uknown>>>>>"; }
    }
    
    ~log_intermediate() noexcept
    {
      try
      {
        //stringstream s;
        //s << std::left << std::setw(30) << function << " " << line
        //  << prefix
        //  << buffer.str()
        //  << " "
        //  << suffix
        //  << flush;
        //s << source << '\t' << line << '\t' << function
        //  << prefix
        //  << buffer.str()
        //  << "\t"
        //  << suffix
        //  << endl
        //  ;
        //UTIL_LOG_EASY_LOGGER_OUT << s.str();
        std::cout << std::left << std::setw(30) << function << " " << line << prefix << buffer.str() << " " << suffix << endl;
        //UTIL_LOG_EASY_LOGGER_OUT << s.str() << UTIL_LOG_EASY_LOGGER_FLUSH;
      }
      catch ( const exception& e )
      { cerr << "\n\n<<<<<\nexception on " << __PRETTY_FUNCTION__ << "\nwhat=" << e.what() << "\n>>>>>\n\n"; }
      catch ( ... )
      { cerr << "\n\n<<<<<\nexception on " << __PRETTY_FUNCTION__ << "\nunknown\n>>>>>\n\n"; }
    }
  };

#endif

}
