package ubic.gemma.script.example

import gemma.gsec.authentication.ManualAuthenticationService
import org.springframework.context.ApplicationContext
import org.springframework.security.core.context.SecurityContextHolder
import ubic.gemma.core.context.EnvironmentProfiles
import ubic.gemma.core.context.SpringContextUtils

class SpringSupport {
    ApplicationContext ctx;

    SpringSupport() {
        ctx = SpringContextUtils.getApplicationContext(new String[]{
                EnvironmentProfiles.PRODUCTION
        });
        authenticateAnonymously()
    }

    SpringSupport(String s1, String s2) {
        ctx = SpringContextUtils.getApplicationContext(new String[]{
                EnvironmentProfiles.PRODUCTION
        });
        authenticate(s1, s2)
    }

    <T> T getBean(Class<T> beanClass) {
        return ctx.getBean(beanClass);
    }

    void authenticate(String s1, String s2) {
        def mas = ctx.getBean(ManualAuthenticationService.class)
        def auth = mas.authenticate(s1, s2)
        def securityContext = SecurityContextHolder.createEmptyContext()
        securityContext.setAuthentication(auth)
        SecurityContextHolder.setContext(securityContext)
    }

    void authenticateAnonymously() {
        def mas = ctx.getBean(ManualAuthenticationService.class)
        def auth = mas.authenticateAnonymously()
        def securityContext = SecurityContextHolder.createEmptyContext()
        securityContext.setAuthentication(auth)
        SecurityContextHolder.setContext(securityContext)
    }
}
