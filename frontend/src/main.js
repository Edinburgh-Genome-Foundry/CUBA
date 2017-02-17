import Vue from 'vue'
import VueRouter from 'vue-router'
import VueResource from 'vue-resource'
import ElementUI from 'element-ui'
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
Vue.use(ElementUI)
Vue.use(widgets)

// Check the users auth status when the app starts
// auth.checkAuth()

const routes = [{
  path: '/home',
  component: Home
}, {
  path: '/about',
  component: About
}, {
  path: '/login',
  component: Login
}
]

scenarios.list.forEach(function (scenario) {
  routes.push({ path: '/' + scenario.infos.path, component: scenario })
})
routes.push({
  path: '*',
  component: Home
})

console.log(routes)

const router = new VueRouter({
  routes
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
