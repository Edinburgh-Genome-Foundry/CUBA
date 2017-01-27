import Vue from 'vue'
import VueRouter from 'vue-router'
import VueResource from 'vue-resource'
import App from './App'
import Home from './pages/Home'
import About from './pages/About'
import Login from './auth/Login'
import scenarios from './pages/scenarios/scenarios'
import ElementUI from 'element-ui'

// import auth from './auth'
// import BootstrapVue from 'bootstrap-vue'

// Globally register bootstrap-vue components

Vue.use(VueResource)
Vue.use(VueRouter)
Vue.use(ElementUI)

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
