<template lang="pug">
.file-uploader
  el-upload(
            :file-list='fileList',
            :on-change='onNew',
            :on-remove='onRemove',
            :auto-upload="false",
            :multiple='multiple',
            action='',
            drag)
    i.el-icon-upload
    .el-upload__text Drop a {{fileDescription}} file here or <em>click to upload</em>
    .el-upload__tip(slot="tip") {{tip}}
  .remove-all(v-if='fileList.length > 1 && multiple')
    el-tooltip(content='Remove all')
      el-button(icon='el-icon-delete',
                @click='removeAll', style='margin-top: 0.6em' circle)
</template>
<script>
export default {
  props: {
    value: {default: () => ([])},
    multiple: {default: true},
    fileDescription: {default: ''},
    tip: {default: ''}
  },
  data () {
    return {
      fileList: [].concat(this.multiple ? this.value : (this.value ? [this.value] : []))
    }
  },
  methods: {
    onNew (file, fileList) {
      if (!this.multiple) {
        this.fileList.pop()
      }
      var self = this
      var reader = new window.FileReader()
      reader.onloadend = function (ev) {
        if (ev.target.readyState === window.FileReader.DONE) {
          file.content = ev.target.result
          file.status = 'success'
          self.fileList.push(file)
          self.$nextTick()
        }
      }
      reader.readAsDataURL(file.raw)
    },
    onRemove (file, fileList) {
      if (this.multiple) {
        var self = this
        var files = [].concat(this.fileList)
        this.removeAll()
        files.map(function (f) {
          if (f.uid !== file.uid) {
            self.fileList.push(f)
          }
        })
      } else {
        this.fileList.pop()
      }
    },
    removeAll () {
      while (this.fileList.length) {
        this.fileList.pop()
      }
    }
  },
  watch: {
    fileList: {
      deep: true,
      handler (val) {
        var returned = this.multiple ? val : (val.length ? val[0] : null)
        this.$emit('input', returned)
      }
    },
    value: {
      deep: true,
      handler (val) {
        // this.fileList = [].concat(this.multiple ? this.value : [this.value])
        if (this.multiple) {
          this.fileList = val
        } else {
          if (val.content) {
            this.fileList = [val]
          } else if (this.fileList.length) {
            this.fileList = val.content ? [val] : []
          }
        }
      }
    }
  }
}
</script>

<style lang='scss' scoped>
.file-uploader {
  display: block;
  text-align: center;
  margin-bottom: 1em;
  // height: auto;
}
</style>
