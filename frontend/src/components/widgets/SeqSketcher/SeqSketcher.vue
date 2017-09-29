<template lang='pug'>
.seq-sketcher
  .file-operations
    //- el-button(icon='document' size='mini' @click='downloadExcel').file-link Download Excel
    el-button.file-link(size='mini' @click='downloadJson') <icon name='file-text-o'></icon> Download JSON
    el-button.file-link(size='mini' @click='showFileDialog=true') <icon name='upload'></icon> Upload Sketches
    el-dialog.part-selector(title="Upload a schema", :visible.sync="showFileDialog", size='small')
      files-uploader(v-model='file', :showSelected='false', :multiple='false')
  textarea.title(v-model='sketchesData.title', type="textarea" rows=1, placeholder='(Enter a title here)')
  p
    input.note(v-model='sketchesData.note', placeholder='(Add some note here)')
  .constructs
    construct(v-for="construct, i in sketchesData.constructs", v-model='sketchesData.constructs[i]',
              @duplicate='duplicateConstruct(i)', @new='newConstruct(i)',
              @delete='deleteConstruct(i)', @moveUp='moveUpConstruct(i)',
              @moveDown='moveDownConstruct(i)',
              :key='construct.id')
</template>

<script>
import construct from './Construct'
import download from 'downloadjs'

export default {
  name: 'seq-sketcher',
  props: {
    value: {
      default: () => ({
        title: '',
        note: '',
        constructs: [
          {
            name: '',
            id: 0,
            parts: [
              {
                category: 'promoter',
                label: 'prom1',
                id: 0
              },
              {
                category: 'CDS',
                label: 'gene with a very very long name',
                id: 1
              },
              {
                category: 'terminator',
                label: 'Term1',
                id: 2
              }
            ]
          }
        ]
      })
    }
  },
  data: function () {
    return {
      sketchesData: this.value,
      constructIdCounter: 0,
      file: {},
      showFileDialog: false
    }
  },
  components: {
    construct
  },
  methods: {
    downloadJson: function () {
      download(JSON.stringify(this.sketchesData, null, ' '), 'design.json')
    },
    newConstruct: function (i) {
      var newConstructs = this.sketchesData.constructs.slice()
      this.constructIdCounter++
      newConstructs.splice(i + 1, 0, {
        name: '',
        note: '',
        id: this.constructIdCounter,
        parts: [
          {
            category: 'user-defined',
            label: '',
            id: 0
          }
        ]
      })
      this.$set(this.sketchesData, 'constructs', newConstructs)
    },
    duplicateConstruct: function (i) {
      this.constructIdCounter++
      var newConstructs = this.sketchesData.constructs.slice()
      var newConstruct = Object.assign({}, JSON.parse(JSON.stringify(newConstructs[i])),
                                       {id: this.constructIdCounter})
      this.constructIdCounter++
      newConstructs.splice(i + 1, 0, newConstruct)
      this.$set(this.sketchesData, 'constructs', newConstructs)
    },
    deleteConstruct: function (i) {
      var newConstructs = this.sketchesData.constructs.slice()
      newConstructs.splice(i - 1, 1)
      this.$set(this.sketchesData, 'constructs', newConstructs)
      console.log('delete', i)
    },
    moveUpConstruct: function (i) {
      if (i === 0) {
        return
      }
      var newConstructs = this.sketchesData.constructs.slice()
      var construct = this.sketchesData.constructs[i]
      newConstructs.splice(i, 1)
      newConstructs.splice(i - 1, 0, construct)
      this.$set(this.sketchesData, 'constructs', newConstructs)
    },
    moveDownConstruct: function (i) {
      if (i === this.sketchesData.constructs.length - 1) {
        return
      }
      this.moveUpConstruct(i + 1)
    }
  },
  watch: {
    'file.content': {
      deep: true,
      handler: function (newval) {
        this.showFileDialog = false
        if (newval) {
          var json = JSON.parse(atob(newval.split(',')[1]))
          this.sketchesData = json
        }
      }
    },
    'sketchesData': {
      deep: true,
      handler: function (newval) {
        this.$emit('input', newval)
      }
    }
  }
}

</script>

<style lang='scss' scoped>
.seq-sketcher {
  width: 100%;
  margin-top: 2em;
  input, textarea {
    border: none !important;
    outline: 0 !important;
    border-color: Transparent;
    font-family: 'Raleway';
    font-size: 1em;
    resize: none;
    &.title {
      font-size: 2em;
      margin-bottom: 0.1em;
      width: 100%;
      text-align: center;

    }
    &.note {
      width: 100%;
    }
  }
  .file-operations {
  text-align: center;
  .file-link {
    display: inline;
    font-size: 1em;
    border: none;
    margin-right: 1em;
  }
  margin-bottom: 2.5em;
}
}
</style>
