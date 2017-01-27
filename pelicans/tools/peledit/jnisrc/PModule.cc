#include <PEL_assertions.hh>
#include "PModule.h"
#include <string.h>

#include <jni.h>
#include <iostream>
#include <stack>
#include <deque>
#include <map>
#include <sstream>
#include <fstream>
#include <iosfwd>

#include <PEL.hh>
#include <PEL_Context.hh>
#include <PEL_Exec.hh>
#include <PEL_Data.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_Expression.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleIterator.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ModulePattern.hh>
#include <PEL_KeywordDataIterator.hh>
#include <PEL_Variable.hh>
#include <PEL_Root.hh>
#include <PEL_Vector.hh>
#include <stringVector.hh>
#include <PEL_Exceptions.hh>

using std::stack ;
using std::string ;


// Internationalization
jclass Class_java_lang_String ;
jmethodID MID_String_getBytes ;
jmethodID MID_String_init ;


void JNU_ThrowByName(JNIEnv *env, const char *name, const char *msg) {
   jclass cls = env->FindClass(name);
   /* if cls is NULL, an exception has already been thrown */
   if (cls != NULL) {
      env->ThrowNew(cls, msg);
   }
   /* free the local ref */
   env->DeleteLocalRef(cls);
}

char* JNU_GetStringNativeChars(JNIEnv* env, jstring jstr) {

   jbyteArray bytes = 0;
   char* result = 0;

   if (env->EnsureLocalCapacity(2) < 0) {
      return 0;  // out of memory error
   }

   bytes = (jbyteArray) env->CallObjectMethod(jstr, MID_String_getBytes);

   //on utilise ExceptionCheck au lieu de ExceptionOccured...
   jboolean exc = env->ExceptionCheck();

   if (!exc) {
      jint len = env->GetArrayLength(bytes);
      result = new char [len+1];
      if (result == 0) {
         JNU_ThrowByName(env, "java/lang/OutOfMemoryError", 0);
         env->DeleteLocalRef(bytes);

         return 0;
      }
      env->GetByteArrayRegion(bytes, 0, len, (jbyte*)result);
      result[len] = 0; // NULL-terminate
   }
   else {
      printf("Exception occured...\n");
   }
   env->DeleteLocalRef(bytes);
   return result;
}

// Native to Unicode caracter translation
jstring JNU_NewStringNative(JNIEnv *env, const char *str)
{
   jstring result = NULL ;
   jbyteArray bytes = 0;
   size_t len;
   if (env->EnsureLocalCapacity(2) < 0) {
      return NULL; /* out of memory error */
   }
   len = strlen(str);
   bytes = env->NewByteArray(len);
   if (bytes != NULL) {
      env->SetByteArrayRegion(bytes, 0, len,
			      reinterpret_cast<jbyte *>(const_cast<char *>(str)));
      result = (jstring)env->NewObject( Class_java_lang_String, MID_String_init, bytes);
      env->DeleteLocalRef( bytes);
   }

   return result;
}

#define GETCLASS(v,n) { v=(jclass)env->NewGlobalRef( env->FindClass(n) )  ;\
                        PEL_ASSERT(v!=NULL) ;}

#define TRACE(level) if (level>=verbosity) std::cerr
// << __FILE__ << ":" << __LINE__ << " "

class Linker {
   public:
      typedef unsigned int Id;
   private:
      static std::map<Id, class IPEL_Module*> id_to_ptr;
      static std::map<class PEL_Module*, Id> ptr_to_id;
      static Id counter;
   public:
      static Id next() {return counter++ ;}
      static void attach(Id id, class IPEL_Module*ptr);
      static void detach(Id id, class IPEL_Module*ptr);
      static class IPEL_Module *get_ptr(Id id) {
         PEL_LABEL( "Linker::IPEL_Module" ) ;
         PEL_CHECK_PRE( id < counter ) ;
         return id_to_ptr[id];
      }
      static Id get_id(class PEL_Module *ptr) {
         PEL_LABEL( "Linker::get_id" ) ;
         Id result = ptr_to_id[ptr] ;

         PEL_CHECK_POST( result < counter ) ;

         return result;
      }
};

class IPEL_Node {
   protected:
      class IPEL_Node *parent;
      IPEL_Node() : parent(0) {}
   public:
      void set_parent(IPEL_Node *m) {parent=m;}
      IPEL_Node *get_parent() const {return parent;}
};

class IPEL_Data : public IPEL_Node {
   protected:
      std::string const& my_keyword ;
      PEL_Data const* my_data ;
   public:
      IPEL_Data(std::string const& key, PEL_Data const* a_data)
            : my_keyword(key), my_data(a_data) {
         PEL_LABEL( "IPEL_Data::IPEL_Data" ) ;
         PEL_CHECK_PRE( a_data!=0 ) ;
         PEL_CHECK_PRE( !key.empty() ) ;
      }
      std::string const& keyword() const {
         PEL_LABEL("IPEL_Data ::keyword") ;
         return my_keyword;
      }
      PEL_Data const* data() const {
         PEL_LABEL("IPEL_Data ::data") ;
         return my_data;
      }
};

class IPEL_Module : public IPEL_Node {
   private:
      Linker::Id id;
   protected:
      std::deque<IPEL_Module *> list_modules;
      std::deque<IPEL_Data *>  list_data;
      std::deque<IPEL_Data *>  list_var;
      IPEL_Module *parent;
      PEL_Module *source;

   public:
      IPEL_Module(PEL_Module *);
      ~IPEL_Module();
      int nb_modules() const {return list_modules.size();}
      int nb_data() const {return list_data.size();}
      int nb_variable() const {return list_var.size();}
      PEL_Module *get_source() const {return source;}
      Linker::Id identifier() {return id;}
      void add_module(IPEL_Module *m) {
         PEL_ASSERT(m->get_parent() == 0);
         list_modules.push_back(m);
         m->set_parent(this);
      }
      void add_data(IPEL_Data *e) {
         PEL_ASSERT(e->get_parent() == 0);
         list_data.push_back(e);
         e->set_parent(this);
      }
      void add_variable(IPEL_Data *e) {
         PEL_ASSERT(e->get_parent() == 0);
         list_var.push_back(e);
         e->set_parent(this);
      }
      IPEL_Module *get_module(int idx) {
         PEL_ASSERT(idx >= 0 && idx < nb_modules());
         return list_modules[idx];
      }
      IPEL_Data  *get_data(int idx) {
         PEL_LABEL( "IPEL_Module::get_data" ) ;

         PEL_ASSERT(idx >= 0 && idx < nb_data());
         return list_data[idx];
      }
      IPEL_Data  *get_variable(int idx) {
         PEL_LABEL( "IPEL_Module::get_variable" ) ;

         PEL_ASSERT(idx >= 0 && idx < nb_variable());
         return list_var[idx];
      }
};

Linker::Id Linker::counter = 0;
std::map<Linker::Id, class IPEL_Module*> Linker::id_to_ptr;
std::map<class PEL_Module*, Linker::Id> Linker::ptr_to_id;
void Linker::attach(Id id, class IPEL_Module*ptr) {
   id_to_ptr[id] = ptr;
   ptr_to_id[ptr->get_source()] = id;
}
void Linker::detach(Id id, class IPEL_Module*ptr) {
   id_to_ptr[id] = 0;
   ptr_to_id[ptr->get_source()] = ~0;
}

IPEL_Module::IPEL_Module(PEL_Module *module)
   : id(Linker::next())
   , source(module) {
   PEL_LABEL( "IPEL_Module::IPEL_Module" ) ;
   PEL_CHECK_PRE( module!=0 ) ;

   Linker::attach(id, this);

   PEL_ModuleIterator *mit = source->create_module_iterator(0);
   for (mit->start() ; mit->is_valid() ; mit->go_next()) {
      add_module(new IPEL_Module(mit->item()));
   }
   mit->destroy();

   PEL_KeywordDataIterator* kwit = source->create_entry_iterator(0);
   for (kwit->start() ; kwit->is_valid() ; kwit->go_next()) {
      add_data(new IPEL_Data(kwit->item()->keyword(),kwit->item()->data()));
   }
   kwit->destroy();

   PEL_Context const* ctx = source->context() ;
   PEL_Module const* father = source->father() ;
   PEL_Context const* father_ctx = ( father!=0 ? father->context() : 0 ) ;


   for (size_t idx=0 ; idx<ctx->nb_variables() ; idx++ ) {
      PEL_Variable const* var = ctx->variable(idx) ;
      if( father_ctx==0 || !father_ctx->has_variable( var ) )
      {
         PEL_Data const* data = ctx->value(var) ;
         add_variable(new IPEL_Data(var->name(),data));
      }
   }
}

IPEL_Module::~IPEL_Module() {
   PEL_LABEL( "IPEL_Module::~IPEL_Module" ) ;

   Linker::detach(id, this);

   for (int idx=0 ; idx < this->nb_variable() ; idx++ ) {
      delete this->get_variable(idx);
   }

   for (int idx=0 ; idx < this->nb_data() ; idx++ ) {
      delete this->get_data(idx);
   }

   for (int idx=0 ; idx < this->nb_modules() ; idx++ ) {
      delete this->get_module(idx);
   }
}



//-----------------------------------------------------------------------------------
// The PData class.
typedef struct {
    jfieldID name;
    jfieldID type;
    jfieldID value;
    jmethodID init;
} DataIdFields;
jclass classDataId = 0;
DataIdFields dataFields;
//    env->SetObjectField(jdata, dataFields.value, p);
//    env->SetObjectField(jdata, dataFields.name, jname);
//    env->SetObjectField(jdata, dataFields.type, jtype);

//-----------------------------------------------------------------------------------
// logging levels used with TRACE
jint verbosity = 0;
jint warning_level = 0;
jint severe_level = 0;
jint config_level = 0;
jint info_level = 0;
jint fine_level = 0;
jclass classLevelId = 0;
jmethodID int_value = 0;

int get_level_value(JNIEnv *env, const char * name) {
   jfieldID fid=env->GetStaticFieldID(classLevelId, name, "Ljava/util/logging/Level;");
   PEL_ASSERT(fid != 0);
   jobject level = env->GetStaticObjectField(classLevelId, fid);
   PEL_ASSERT(level != 0);
   return env->CallIntMethod(level, int_value);
}
//-----------------------------------------------------------------------------------
// Initializes the whole package : must be called only once.
JNIEXPORT void JNICALL Java_data_jni_PModule_initialize(JNIEnv *env, jclass, jint v) {
   PEL_LABEL( "JNICALL::initialize" ) ;

   GETCLASS( Class_java_lang_String , "java/lang/String" );
   MID_String_getBytes = env->GetMethodID(Class_java_lang_String, "getBytes", "()[B");
   MID_String_init = env->GetMethodID(Class_java_lang_String, "<init>", "([B)V");

   PEL_ASSERT(classLevelId == 0);
   GETCLASS(classLevelId, "java/util/logging/Level");
   int_value = env->GetMethodID(classLevelId, "intValue", "()I");
   PEL_ASSERT(int_value != 0);

   verbosity = v;
   severe_level = get_level_value(env, "SEVERE");
   warning_level = get_level_value(env, "WARNING");
   config_level = get_level_value(env, "CONFIG");
   info_level = get_level_value(env, "INFO");
   fine_level = get_level_value(env, "FINE");
   TRACE(fine_level) << " PEL_Exec::initialize start verbosity=" << verbosity << std::endl;

   PEL_ASSERT(classDataId == 0);
   GETCLASS(classDataId, "data/jni/PData");

   dataFields.name  = env->GetFieldID(classDataId, "name", "Ljava/lang/String;");
   PEL_ASSERT( dataFields.name!=0 ) ;
   dataFields.type  = env->GetFieldID(classDataId, "type", "Ljava/lang/String;");
   PEL_ASSERT( dataFields.type!=0 ) ;
   dataFields.value = env->GetFieldID(classDataId, "value", "Ljava/lang/String;");
   PEL_ASSERT( dataFields.value!=0 ) ;
   dataFields.init  = env->GetMethodID(classDataId, "<init>", "()V");
   PEL_ASSERT( dataFields.init!=0 ) ;

   // PEL initialize
   stringVector args(0) ;
   char *argv[] = {
      // process name
      const_cast<char*>("JNI")
      // no signal handling incompatible with Java
      , const_cast<char*>("-no_signal_handling")
      // no external API (MPI nor PETSc)
      , const_cast<char*>("-no_external_API")
   };
   int ok = PEL_Exec::initialize(3, argv, args) ;
   TRACE(fine_level) << " PEL_Exec::initialize returns " << ok << std::endl;

}

//-----------------------------------------------------------------------------------
// Loads a Module from file ans returns the identifier.
JNIEXPORT jint JNICALL
Java_data_jni_PModule_createFromFile__Ljava_lang_String_2(JNIEnv *env, jclass, jstring name) {
   PEL_LABEL( "JNICALL::createFromFile" ) ;
   const char *filename = JNU_GetStringNativeChars(env, name);
   PEL_ASSERT( filename!=0 ) ;

   //TRACE << "(" << filename << ")" << std::endl;
   int id = -1;
   try {
      PEL_Module *a_mod = PEL_Module::create( 0, "root", filename ) ;
      //std::cout << "Fichier lu \n" << std::endl ;
      //std::cout.flush() ;

      IPEL_Module *ipel = new IPEL_Module(a_mod);
      id = ipel->identifier();

   } catch (PEL_Exceptions::Error) {
      TRACE(severe_level) << " Exception raised" << std::endl;
   }

   //TRACE << filename << " read id=" << id << std::endl;
   PEL_CHECK( !env->ExceptionCheck() ) ;

   return id;
}
//-----------------------------------------------------------------------------------
// Loads a Module from file at position pos and returns the identifier.
JNIEXPORT jint JNICALL
Java_data_jni_PModule_createFromFile__Ljava_lang_String_2J(JNIEnv *env, jclass, jstring name, jlong pos) {
   PEL_LABEL( "JNICALL::createFromFile&position" ) ;
   const char *filename = JNU_GetStringNativeChars(env, name);
   PEL_ASSERT( filename!=0 ) ;

   TRACE(fine_level) << "(" << filename << ", " << pos << ")" << std::endl;
   int id = -1;
   try {
      std::ifstream file( filename ) ;
      file.seekg(pos);

      PEL_Module *a_mod = PEL_Module::create( 0, "root", file ) ;
      std::cout << "Fichier lu \n" << std::endl ;
      std::cout.flush() ;

      IPEL_Module *ipel = new IPEL_Module(a_mod);
      id = ipel->identifier();

   } catch (PEL_Exceptions::Error) {
      TRACE(severe_level) << " Exception raised" << std::endl;
   }

   TRACE(fine_level) << filename << ", " << pos << " read id=" << id << std::endl;
   delete [] filename;
   PEL_CHECK( !env->ExceptionCheck() ) ;

   return id;
}

//-----------------------------------------------------------------------------------
// Loads a Module from a string ans returns the identifier.
JNIEXPORT jint JNICALL Java_data_jni_PModule_createFromString(JNIEnv *env, jclass, jstring str) {
   PEL_LABEL( "JNICALL::createFromString" ) ;
   const char *string = JNU_GetStringNativeChars(env, str);
   PEL_ASSERT( string!=0 ) ;

   std::istringstream is( string ) ;
   int id = -1;
   try {
      PEL_Module *a_mod = PEL_Module::create( 0, "root", is ) ;
      IPEL_Module *ipel = new IPEL_Module(a_mod);
      id = ipel->identifier();

   } catch (PEL_Exceptions::Error) {
      TRACE(severe_level) << " Exception raised" << std::endl;
   }

   TRACE(fine_level) << "read id=" << id << std::endl;
   delete [] string ;
   PEL_CHECK( !env->ExceptionCheck() ) ;

   return id;
}
//-----------------------------------------------------------------------------------
// returns the new string containing the name of the module pointed by id
JNIEXPORT jstring JNICALL
Java_data_jni_PModule_name(JNIEnv *env, jclass, jint id) {
    PEL_LABEL( "Java_data_jni_PModule_name" ) ;
    IPEL_Module * ipel = Linker::get_ptr(id);
    jstring result = JNU_NewStringNative( env, ipel->get_source()->name().c_str() );
    PEL_CHECK_POST( !env->ExceptionCheck() ) ;

    return result;
}
//-----------------------------------------------------------------------------------
// deletes the module pointed by id
JNIEXPORT void JNICALL
Java_data_jni_PModule_destroy(JNIEnv *env, jclass, jint id) {
    PEL_LABEL( "Java_data_jni_PModule_name" ) ;
    IPEL_Module * ipel = Linker::get_ptr(id);
    PEL_ASSERT(ipel != 0);
    PEL_ASSERT(ipel->get_parent() == 0);

    ipel->get_source()->destroy();
    delete ipel;

    PEL_CHECK_POST( !env->ExceptionCheck() ) ;
}


//-----------------------------------------------------------------------------------
// returns the number of data in the module pointed by id
JNIEXPORT jint JNICALL
Java_data_jni_PModule_nbData(JNIEnv *, jclass, jint id) {
    PEL_LABEL( "Java_data_jni_PModule_nbData" ) ;
    IPEL_Module * ipel = Linker::get_ptr(id);
    PEL_ASSERT(ipel != 0);
    return ipel->nb_data();
}


//-----------------------------------------------------------------------------------
// returns the number of modules in the module pointed by id
JNIEXPORT jint JNICALL
Java_data_jni_PModule_nbModules(JNIEnv *, jclass, jint id) {
   PEL_LABEL( "Java_data_jni_PModule_nbModules" );
   IPEL_Module * ipel = Linker::get_ptr(id);
   PEL_ASSERT(ipel != 0);
   return ipel->nb_modules();
}

//-----------------------------------------------------------------------------------
// returns the identifier of the ieth submodule in the module pointed by id
JNIEXPORT jint JNICALL
Java_data_jni_PModule_getModule__II(JNIEnv *, jclass, jint id, jint ieth) {
   PEL_LABEL( "JNICALL::getModule" ) ;
   jint result = -1 ;

   IPEL_Module * ipel = Linker::get_ptr(id);
   if( ipel!=0 ) result=ipel->get_module(ieth)->identifier() ;

   return result;
}

//-----------------------------------------------------------------------------------
// returns the identifier of the submodule named 'name' in the module pointed by id
JNIEXPORT jint JNICALL
Java_data_jni_PModule_getModule__ILjava_lang_String_2(JNIEnv *env, jclass, jint id, jstring name) {
   PEL_LABEL( "JNICALL::getModule" ) ;
   IPEL_Module * ipel = Linker::get_ptr(id);
   PEL_ASSERT(ipel != 0);

   const char *str = JNU_GetStringNativeChars(env, name);
   PEL_ASSERT( str!=0 ) ;
   int result_id = -1;
   if( ipel->get_source()->has_module(str) )
   {
      try {
         PEL_Module *sub = ipel->get_source()->module(str);
         result_id = Linker::get_id(sub);
         PEL_ASSERT(result_id != 0);
      } catch (PEL_Exceptions::Error) {
         TRACE(severe_level) << " Exception raised" << std::endl;
      }
   }

   delete [] str;
   PEL_CHECK( !env->ExceptionCheck() ) ;

   return result_id;
}

jstring
to_string(JNIEnv *env, const PEL_Data *data) {
   PEL_LABEL( "to_string" ) ;

   std::ostringstream st;
   st.precision( 15 ) ;

   data->print( st, 0 ) ;
   jstring result = JNU_NewStringNative( env, st.str().c_str() );

   PEL_CHECK_POST( !env->ExceptionCheck() ) ;
   return result ;

}

//-----------------------------------------------------------------------------------
// sets the data object to the value of the ieth data in the module pointed by id
JNIEXPORT jint JNICALL
Java_data_jni_PModule_getData__IILdata_jni_PData_2(JNIEnv *env, jclass, jint id, jint ieth, jobject jdata) {
   PEL_LABEL( "JNICALL::getData" ) ;
   IPEL_Module * ipel = Linker::get_ptr(id);
   PEL_ASSERT(ipel != 0);
   IPEL_Data * idata =  ipel->get_data(ieth);
   PEL_ASSERT(idata != 0);
   std::string typname = PEL_Data::type_name( idata->data()->data_type() ) ;
   env->SetObjectField(jdata, dataFields.name, JNU_NewStringNative(env, idata->keyword().c_str()));
   env->SetObjectField(jdata, dataFields.type, JNU_NewStringNative(env, typname.c_str()));
   env->SetObjectField(jdata, dataFields.value, to_string(env, idata->data()));
   PEL_CHECK( !env->ExceptionCheck() ) ;
   return 1;
}

//-----------------------------------------------------------------------------------
// sets the data object to the value of the data named 'name' in the module pointed by id
JNIEXPORT jint JNICALL
Java_data_jni_PModule_getData__ILjava_lang_String_2Ldata_jni_PData_2(JNIEnv *env, jclass, jint id, jstring name, jobject jdata) {
   PEL_LABEL( "JNICALL::getData" ) ;
   IPEL_Module * ipel = Linker::get_ptr(id);
   PEL_ASSERT(ipel != 0);
   const char *str = JNU_GetStringNativeChars(env, name);
   PEL_ASSERT( str!=0 ) ;

   int result = 0;
   try {
      if (ipel->get_source()->has_entry(str)) {
	 const PEL_Data *data = ipel->get_source()->data_of_entry(str);
	 result = 1;
	 env->SetObjectField(jdata, dataFields.name, name);
	 std::string typname = PEL_Data::type_name( data->data_type() ) ;
	 env->SetObjectField(jdata, dataFields.type, JNU_NewStringNative(env, typname.c_str()));
	 env->SetObjectField(jdata, dataFields.value, to_string(env, data));
	 PEL_CHECK( !env->ExceptionCheck() ) ;
      }
   } catch (PEL_Exceptions::Error) {
      TRACE(severe_level) << " Exception raised" << std::endl;
   }
   delete []( str);
   PEL_CHECK( !env->ExceptionCheck() ) ;
   return result;
}

JNIEXPORT jobject JNICALL
Java_data_jni_PModule_ngetData(JNIEnv *env, jclass jcl, jint id, jint ieth) {
   PEL_LABEL( "JNICALL::ngetData" ) ;
   PEL_CHECK_PRE(env != 0);
   PEL_CHECK_PRE(id != 0);
   PEL_CHECK_PRE( ieth <
                  Java_data_jni_PModule_nbData(env, jcl, id ) ) ;
   IPEL_Module * ipel = Linker::get_ptr(id);
   PEL_ASSERT(ipel != 0);
   IPEL_Data * idata = ipel->get_data(ieth);
   PEL_ASSERT(idata != 0);

   jobject jdata = env->NewObject( classDataId,  dataFields.init ) ;
   PEL_CHECK( !env->ExceptionCheck() ) ;


   std::string typname = PEL_Data::type_name( idata->data()->data_type() ) ;
   env->SetObjectField(jdata, dataFields.name, JNU_NewStringNative(env, idata->keyword().c_str()));
   PEL_CHECK( !env->ExceptionCheck() ) ;
   env->SetObjectField(jdata, dataFields.type, JNU_NewStringNative(env, typname.c_str()));
   PEL_CHECK( !env->ExceptionCheck() ) ;
   env->SetObjectField(jdata, dataFields.value, to_string(env, idata->data()));
   PEL_CHECK( !env->ExceptionCheck() ) ;

   return jdata;
}


//-----------------------------------------------------------------------------------
// returns the number of variables in the module pointed by id
JNIEXPORT jint JNICALL
Java_data_jni_PModule_nbVariable(JNIEnv *, jclass, jint id) {
    PEL_LABEL( "Java_data_jni_PModule_nbVariable" ) ;
    IPEL_Module * ipel = Linker::get_ptr(id);
    PEL_ASSERT(ipel != 0);
    return ipel->nb_variable();
}

//-----------------------------------------------------------------------------------
// returns the ieth variable
JNIEXPORT jobject JNICALL
Java_data_jni_PModule_ngetVariable(JNIEnv *env, jclass jcl, jint id, jint ieth) {
   PEL_LABEL( "JNICALL::ngetVariable" ) ;
   PEL_CHECK_PRE(env != 0);
   PEL_CHECK_PRE(id != 0);
   PEL_CHECK_PRE( ieth <
                  Java_data_jni_PModule_nbVariable(env, jcl, id ) ) ;
   IPEL_Module * ipel = Linker::get_ptr(id);
   PEL_ASSERT(ipel != 0);
   IPEL_Data * idata = ipel->get_variable(ieth);
   PEL_ASSERT(idata != 0);

   jobject jdata = env->NewObject( classDataId,  dataFields.init ) ;
   PEL_CHECK( !env->ExceptionCheck() ) ;


   std::string typname = PEL_Data::type_name( idata->data()->data_type() ) ;
   env->SetObjectField(jdata, dataFields.name, JNU_NewStringNative(env, idata->keyword().c_str()));
   PEL_CHECK( !env->ExceptionCheck() ) ;
   env->SetObjectField(jdata, dataFields.type, JNU_NewStringNative(env, typname.c_str()));
   PEL_CHECK( !env->ExceptionCheck() ) ;
   env->SetObjectField(jdata, dataFields.value, to_string(env, idata->data()));
   PEL_CHECK( !env->ExceptionCheck() ) ;

   return jdata;
}

/*
 * Class:     jni_PModule
 * Method:    evaluate
 * Signature: (Ljava/lang/String;Ljava/lang/String;I)Z
 */
JNIEXPORT jstring JNICALL Java_data_jni_PModule_evaluate_1entry(JNIEnv *env, jclass, jstring module, jstring key)
{
   PEL_LABEL( "PModule::evaluate" ) ;
   const char *str = JNU_GetStringNativeChars(env, module);
   const char *keystr = JNU_GetStringNativeChars(env, key);

   std::istringstream is( str ) ;
   jstring result = 0 ;
   PEL_Module *a_mod = 0;
   try {
      a_mod = PEL_Module::create( 0, "root", is ) ;
      PEL_ModuleIterator* it = a_mod->create_module_iterator(a_mod) ;
      it->start() ;
      PEL_Module *tested_mod = it->item();
      bool failed = false ;
      std::string eval_by_mod = tested_mod->data_as_string(keystr,0,failed) ;
      if( !failed )
      {
         result = JNU_NewStringNative( env, eval_by_mod.c_str() ) ;
      }
      else
      {
         TRACE(config_level) << " evaluate_entry : can't evaluate expression : "
                             <<keystr << std::endl << " in module :" << std::endl << str << std::endl;
      }

   } catch (PEL_Exceptions::Error) {
      TRACE(severe_level) << " evaluate_entry : Exception raised in JNI evaluate function when evaluating expression : " <<
      keystr << std::endl << " in module :" << std::endl << str << std::endl;
   }

   if (a_mod != 0) a_mod->destroy() ;
   delete [] keystr;
   delete [] str;
   PEL_CHECK( !env->ExceptionCheck() ) ;

   return result ;
}


JNIEXPORT jstring JNICALL Java_data_jni_PModule_evaluate_1expression(JNIEnv *env, jclass, jstring module, jstring expr)
{
   PEL_LABEL( "PModule::evaluate" ) ;
   const char *str = JNU_GetStringNativeChars(env, module);
   const char *keystr = JNU_GetStringNativeChars(env, expr);

   std::istringstream is( str ) ;
   jstring result = 0 ;
   PEL_Module *a_mod = 0;
   try {
      a_mod = PEL_Module::create( 0, "root", is ) ;
      PEL_ModuleIterator* it = a_mod->create_module_iterator(a_mod) ;
      it->start() ;
      PEL_Module *tested_mod = it->item();

      PEL_Data const* val = tested_mod->create_evaluation(a_mod,keystr,0) ;
      if( val != 0 )
      {
         PEL_Data const* simp = val->create_simplification(a_mod) ;
         result = to_string( env, simp ) ;
      }
   } catch (PEL_Exceptions::Error) {
      TRACE(severe_level) << " evaluate_expression : Exception raised in JNI evaluate function when evaluating expression : " <<
      keystr << std::endl << " in module :" << std::endl << str << std::endl;
   }

   if (a_mod != 0) a_mod->destroy() ;
   delete [] keystr;
   delete [] str;
   PEL_CHECK( !env->ExceptionCheck() ) ;

   return result ;
}

/*
 * Class:     data_jni_PModule
 * Method:    getRegisteredExpressions
 * Signature: ()[Ljava/lang/String;
 */
JNIEXPORT jobjectArray JNICALL Java_data_jni_PModule_getRegisteredExpressions(JNIEnv *env , jclass) {
	PEL_LABEL( "PModule::getRegisteredExpressions" ) ;
	PEL_CHECK_PRE(env != 0);

  	stringVector exprs = PEL_Expression::registered_expressions();
  	exprs.sort();

  	jclass stringClass;
  	GETCLASS(stringClass, "java/lang/String");
  	jobjectArray result = env->NewObjectArray(exprs.size(), stringClass, 0);
  	PEL_ASSERT(result != 0);

  	for(unsigned int i = 0; i < exprs.size(); i++) {
  		std::string ieth = PEL_Expression::usage_of(exprs(i));
  		env->SetObjectArrayElement(result, i , JNU_NewStringNative(env, ieth.c_str()));
    	PEL_CHECK( !env->ExceptionCheck() ) ;
  	}

    return result ;
}

