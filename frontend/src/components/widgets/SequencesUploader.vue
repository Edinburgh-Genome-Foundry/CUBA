<template lang="pug">
.sequences-uploader(:class='{withfiles: (files.length != 0)}')
  filesuploader(v-model='files', :help='help', :text='text',
                :multiple='true',
                tip="Accepted formats: FASTA, Genbank, Snapgene")
  p.num-files(v-if='files.length') {{files.length}} {{files.length > 1 ? 'files' : 'file'}} selected
  .sequences-list
    el-row(v-for='(file, index) in filesWithLinearities', :key='index', :gutter='20')
      el-col(:span='6')
        el-select(v-model='circularities[index]', size='mini')
          el-option(label='circular', :value='true')
          el-option(label='linear', :value='false')
      el-col(:span='18') {{file.name}}
</template>

<script>
import filesuploader from './FilesUploader'

export default {
  props: {
    text: {default: 'Drop files here or click to select'},
    help: {default: 'Accepted formats: genbank, fasta'},
    multiple: {default: true},
    filter: {default: () => () => true},
    displaySelected: {default: false}
  },
  data () {
    return {
      files: [],
      circularities: [],
      hovering: false
    }
  },
  computed: {
    filesWithLinearities () {
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
    files (val) {
      var self = this
      self.circularities = []
      val.forEach(function (f) {
        self.circularities.push(true)
      })
    },
    circularities (val) {
      console.log(val)
    },
    filesWithLinearities (val) {
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
