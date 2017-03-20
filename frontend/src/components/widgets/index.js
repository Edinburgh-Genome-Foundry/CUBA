export default function install (Vue) {
  Vue.component('backend-querier', require('./BackendQuerier'))
  Vue.component('learn-more', require('./LearnMore'))
  Vue.component('download-button', require('./DownloadButton'))
  Vue.component('helper', require('./Helper'))
}
