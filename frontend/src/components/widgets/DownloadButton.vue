<template lang='pug'>
el-button.download.center(type='primary' @click='downloadContent' icon='document') {{text}}
</template>

<script>
var download = require('downloadjs')
export default {
  props: {
    text: {
      default: 'Download'
    },
    filedata: {}
  },
  methods: {
    downloadContent () {
      // console.log(this.filedata.mimetype, download)
      download(this.filedata.data, this.filedata.name, this.filedata.mimetype)
    }
  }
}
</script>

<style scoped>
</style>
