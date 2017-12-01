import Vue from 'vue'
import VueRouter from 'vue-router'
import VueResource from 'vue-resource'
import Icon from 'vue-awesome/components/Icon'
import 'vue-awesome/icons'
import ElementUI from 'element-ui'
import 'element-ui/lib/theme-chalk/index.css'
import locale from 'element-ui/lib/locale/lang/en'
import App from './App'
import Home from './components/pages/Home'
import About from './components/pages/About'
import Login from './auth/Login'
import scenarios from './components/scenarios/scenarios'
import widgets from './components/widgets'

// import auth from './auth'
// import BootstrapVue from 'bootstrap-vue'

// Globally register bootstrap-vue components

Vue.use(VueResource)
Vue.use(VueRouter)
Vue.use(widgets)
Vue.component('icon', Icon)
Vue.use(ElementUI, { locale })

// Check the users auth status when the app starts
// auth.checkAuth()

const routes = [{
  path: '/home',
  component: Home,
  meta: {
    title: 'EGF CUBA',
    description: 'The Collection of Useful SynBio apps'
  }
}, {
  path: '/about',
  component: About,
  meta: {title: 'About EGF\'s CUBA'}
}, {
  path: '/login',
  component: Login,
  meta: {title: 'CUBA - Login'}
}
]

scenarios.list.forEach(function (category) {
  category.scenarios.forEach(function (scenario) {
    routes.push({
      path: '/' + scenario.infos.path,
      component: scenario,
      meta: {
        title: scenario.infos.title + '- CUBA',
        description: 'From the EGF\'s Collection of useful (Syn)Bio Apps}'
      }
    })
  })
})

routes.push({
  path: '*',
  component: Home,
  meta: {
    title: 'EGF CUBA',
    description: 'The Collection of Useful SynBio apps'
  }
})

const router = new VueRouter({
  routes,
  mode: 'history'
})
router.beforeEach((to, from, next) => {
  document.title = to.meta.title
  next()
})

/* eslint-disable no-new */
new Vue({
  router,
  el: '#app',
  template: '<App/>',
  components: {
    App
  }
})
