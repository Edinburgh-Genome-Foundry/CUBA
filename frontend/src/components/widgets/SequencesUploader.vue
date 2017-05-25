<template lang="pug">
.sequences-uploader(:class='{withfiles: (files.length != 0)}')
  filesuploader(v-model='files', :help='help', :text='text', :multiple='multiple')
  p.num-files(v-if='files.length') {{files.length}} {{files.length > 1 ? 'files' : 'file'}} selected
  .sequences-list
    el-row(v-for='(file, index) in filesWithLinearities', :key='index', :gutter='40')
      el-col(:span='5')
        el-switch(v-model='circularities[index]', :width='90',
                  on-text='circular', off-text='linear')
      el-col(:span='18') {{file.name}}
</template>

<script>
import filesuploader from './FilesUploader'

export default {
  props: {
    text: {default: 'Drop files here or click to select'},
    help: {default: 'Accepted formats: genbank, fasta'},
    multiple: {default: true},
    filter: {default: () => () => true}
  },
  data: function () {
    return {
      files: [],
      circularities: [],
      hovering: false
    }
  },
  computed: {
    filesWithLinearities: function () {
      var result = []
      var self = this
      this.files.forEach(function (f, i) {
        f.circularity = self.circularities[i]
        result.push(f)
      })
      return result
    }
  },
  watch: {
    files: function (val) {
      var self = this
      self.circularities = []
      val.forEach(function (f) {
        self.circularities.push(true)
      })
    },
    circularities: function (val) {
      console.log(val)
    },
    filesWithLinearities: function (val) {
      this.$emit('input', this.multiple ? val : val[0])
    }
  },
  components: {
    filesuploader
  }
}
</script>

<style scoped>
.sequences-list{
  max-height: 413px;
  overflow-y: auto;
  overflow-x: hidden;
  margin-top: 30px;
  margin-left: 10%;
}
.el-switch {
  margin-bottom: 10px;
}
.el-switch__label.el-switch__label--right {
  width: 100% !important;
}

p.num-files {
  font-size: 14px;
  margin-top: 5px
}

.sequences-uploader {
  padding:10px;
}
.withfiles {
  background-color: rgba(34, 157, 249, 0.03);
  border: 0.5px solid #eeeeee;
  border-radius: 10px;
}
</style>
