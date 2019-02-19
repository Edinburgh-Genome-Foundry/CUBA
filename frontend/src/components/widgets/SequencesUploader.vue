<template lang="pug">
.sequences-uploader(:class='{withfiles: (files.length != 0)}')
  filesuploader(v-model='files', :help='help', :text='text',
                :multiple='true',
                tip="Accepted formats: FASTA, Genbank, Snapgene")
  p.num-files(v-if='files.length') {{files.length}} {{files.length > 1 ? 'files' : 'file'}} selected
  .sequences-list
    el-row(v-for='(file, index) in files', :key='index', :gutter='20')
      el-col(:span='6')
        el-select(v-model='files[index].circularity', size='mini')
          el-option(label='circular', :value='true')
          el-option(label='linear', :value='false')
      el-col(:span='18') {{file.name}}
</template>

<script>
import filesuploader from './FilesUploader'

export default {
  props: {
    value: {default: () => ([])},
    text: {default: 'Drop files here or click to select'},
    help: {default: 'Accepted formats: genbank, fasta'},
    multiple: {default: true},
    filter: {default: () => () => true},
    displaySelected: {default: false}
  },
  data () {
    return {
      files: this.value,
      hovering: false
    }
  },
  watch: {
    value: {
      deep: true,
      handler (val) {
        val.map(function (file) {
          if (file.circularity !== true) {
            file.circularity = false
          }
        })
        this.files = val
      }
    },
    files: {
      deep: true,
      handler (val) {
        this.$emit('input', this.multiple ? val : val[0])
      }
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
