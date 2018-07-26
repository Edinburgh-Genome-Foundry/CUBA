<template lang="pug">
.file-example
  el-dialog(title="", :visible.sync="dialogVisible" width="80%")
    center
      img.big-image(v-if='imgSrc', :src='imgSrc')

  el-card
    center
      img.small-image(v-if='imgSrc', :src='imgSrc' @click="dialogVisible = true")
      h3
        .filename {{dataFilename}}
      el-row(:gutter='60')
        span.use-file(v-if='useFileButton')
          el-tooltip(content='use this file')
            el-button(icon='el-icon-plus' @click='loadFile()' circle)
        span &nbsp;
        el-tooltip(content='download this file')
          a(:href='fileHref')
            el-button(icon='el-icon-download' circle)
    p(v-if='description') {{description}}
    slot
  //- hr
</template>

<script>
import uuidv1 from 'uuid/v1'

export default {
  name: 'file-example',
  props: {
    fileHref: {default: () => ''},
    imgSrc: {default: () => null},
    description: {default: () => null},
    filename: {default: () => null},
    useFileButton: {default: true}
  },
  data () {
    return {
      defaultImgSrc: null,
      dataFilename: this.filename,
      dialogVisible: false
    }
  },
  methods: {
    loadFile () {
      console.log('loadfile')
      var self = this
      this.$http.get(this.fileHref, {responseType: 'blob'}).then(function (response) {
        console.log(response)
        var reader = new window.FileReader()
        reader.onloadend = function (ev) {
          if (ev.target.readyState === window.FileReader.DONE) {
            self.$emit('input', {
              name: response.url.split('/').slice(-1)[0],
              content: ev.target.result,
              uid: uuidv1(),
              status: 'success'
            })
          }
        }
        reader.readAsDataURL(response.bodyBlob)
      })
    }
  }
}
</script>

<style lang='scss' scoped>

.file-example {

  hr {
    margin-bottom: 1em;
    border: 0;
    height: 1px;
    background: #333;
    background-image: linear-gradient(to right, #ccc, #333, #ccc);
  }
  img.small-image {
    max-width: 100%;
    max-height: 9em;
    cursor: pointer;
  }
  img.big-image {
    max-width: 100%;
  }

  .text-column {
    text-align: left;
    h3 {
      text-align: center;
      .fa-icon {
        margin-bottom: -0.2em;
      }
    }
    a {
      &:hover {

      }

      color: black;
      text-decoration: none;
      .filename {
        display: inline-block;
        margin-left: 0.5em;
      }
    }
  }
}
</style>
