<template lang='pug'>
a.btn.btn-lg.btn-default.center-block.download-button(type='button' @click='downloadContent')
  i.el-icon-document
  span Download Report
</template>

<script>
var download = require('downloadjs')
export default {
  props: {
    button_text: {
      default: 'Submit'
    },
    filedata: {}
  },
  methods: {
    downloadContent: function () {
      // console.log(this.filedata.mimetype, download)
      download(this.filedata.data, this.filedata.name, this.filedata.mimetype)
    }
  }
}
</script>

<style scoped>
.download-button {
  margin-top: 20px;
  max-width: 230px;
}
.el-icon-document {
  display: inline;
  margin-right: 10px
}
</style>
