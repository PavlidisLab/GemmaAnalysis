package ubic.gemma.script.example

import gemma.gsec.authentication.ManualAuthenticationService
import org.springframework.beans.BeansException
import org.springframework.beans.factory.BeanFactory
import org.springframework.beans.factory.NoSuchBeanDefinitionException
import org.springframework.context.ApplicationContext
import org.springframework.context.ConfigurableApplicationContext
import org.springframework.security.core.context.SecurityContextHolder
import ubic.gemma.core.context.EnvironmentProfiles
import ubic.gemma.core.context.SpringContextUtils

class SpringSupport implements BeanFactory {

    private ApplicationContext ctx;

    SpringSupport() {
        ctx = SpringContextUtils.getApplicationContext(new String[]{
                EnvironmentProfiles.PRODUCTION
        });
        addShutdownHook {
            shutdown()
        }
        authenticateAnonymously()
    }

    SpringSupport(String username, String password) {
        ctx = SpringContextUtils.getApplicationContext(new String[]{
                EnvironmentProfiles.PRODUCTION
        });
        addShutdownHook {
            shutdown()
        }
        authenticate(username, password)
    }

    @Override
    Object getBean(String name) throws BeansException {
        return ctx.getBean(name)
    }

    @Override
    <T> T getBean(String name, Class<T> beanClass) {
        return ctx.getBean(name, beanClass)
    }

    @Override
    <T> T getBean(Class<T> beanClass) {
        return ctx.getBean(beanClass)
    }

    @Override
    Object getBean(String name, Object... args) throws BeansException {
        return ctx.getBean(name, args)
    }

    @Override
    boolean containsBean(String name) {
        return ctx.containsBean(name)
    }

    @Override
    boolean isSingleton(String name) throws NoSuchBeanDefinitionException {
        return ctx.isSingleton(name);
    }

    @Override
    boolean isPrototype(String name) throws NoSuchBeanDefinitionException {
        return ctx.isPrototype(name);
    }

    @Override
    boolean isTypeMatch(String name, Class<?> targetType) throws NoSuchBeanDefinitionException {
        return ctx.isTypeMatch(name, targetType)
    }

    @Override
    Class<?> getType(String name) throws NoSuchBeanDefinitionException {
        return ctx.getType(name)
    }

    @Override
    String[] getAliases(String name) {
        return ctx.getAliases(name)
    }

    void authenticate(String username, String password) {
        def mas = ctx.getBean(ManualAuthenticationService.class)
        def auth = mas.authenticate(username, password)
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

    /**
     * Explicitly shutdown the Spring context.
     *
     * This is not necessary since a shutdown hook is already registered when the context is created by this class.
     */
    void shutdown() {
        if (ctx instanceof ConfigurableApplicationContext) {
            ctx.close()
        }
    }
}
