#ifndef COND_XMLAUTHENTITACTIONSERVICE_H
#define COND_XMLAUTHENTITACTIONSERVICE_H

#include "CondCore/CondDB/interface/CredentialStore.h"
#include "CondCore/CondDB/src/IDbAuthentication.h"
//
#include "RelationalAccess/IAuthenticationService.h"
#include "CoralKernel/Service.h"
#include "CoralKernel/Property.h"
//
#include <map>
#include <set>
#include <string>

namespace coral {

  class AuthenticationCredentials;
  //class IAuthenticationCredentials;
}  // namespace coral

namespace cond {

  namespace RelationalAuthenticationService {

    /**
     */
    class RelationalAuthenticationService : public coral::Service,
                                            virtual public coral::IAuthenticationService,
                                            virtual public persistency::IDbAuthentication {
    public:
      /// Standard Constructor
      explicit RelationalAuthenticationService(const std::string& name);

      /// Standard Destructor
      ~RelationalAuthenticationService() override;

    public:
      /// Sets the input file name
      void setAuthenticationPath(const std::string& inputPath);

      /**
       * Returns a reference to the credentials object for a given connection string.
       * If the connection string is not known to the service an UnknownConnectionException is thrown.
       */
      const coral::IAuthenticationCredentials& credentials(const std::string& connectionString) const override;

      /**
       * Returns a reference to the credentials object for a given connection string.
       * If the connection string is not known to the service an UnknownConnectionException is thrown.
       * If the role is not known to the service an UnknownRoleException is thrown.
       */
      const coral::IAuthenticationCredentials& credentials(const std::string& connectionString,
                                                           const std::string& role) const override;

      std::string principalName() override;

    private:
      /// The input file with the data
      std::string m_authenticationPath;

      /// The service providing the authentication data
      mutable CredentialStore m_db;

      mutable coral_bridge::AuthenticationCredentialSet m_cache;

      coral::Property::CallbackID m_callbackID;
    };

  }  // namespace RelationalAuthenticationService

}  // namespace cond

#endif
