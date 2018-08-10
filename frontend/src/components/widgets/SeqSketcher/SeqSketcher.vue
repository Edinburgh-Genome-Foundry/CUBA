<template lang='pug'>
.seq-sketcher
  .file-operations
    //- el-button(icon='document' size='mini' @click='downloadExcel').file-link Download Excel
    el-button-group
      el-button(icon='el-icon-upload2' size='mini' @click='showFileDialog=true') Upload sketches
      el-button(icon='el-icon-download' size='mini' @click='downloadJson') Download as JSON
      el-button(icon='el-icon-download' size='mini' @click='downloadXlsx') Download as Xlsx
    el-dialog.part-selector(title="Upload a schema", :visible.sync="showFileDialog", size='small')
      collapsible(title='Examples')
        file-example(filename='example_construct_sketches.xlsx',
                     @input='function (e) {file = e}',
                     fileHref='/static/file_examples/sketch_constructs/example_construct_sketches.xlsx',
                     imgSrc='/static/file_examples/sketch_constructs/example_construct_sketches.png')
          p.
            Collection of constructs assembled using random parts from the EMMA standard.
            Unzip the file and drag the genbank files into the file upload area.
      files-uploader(v-model='file', :showSelected='false', :multiple='false')
  textarea.title(v-model='sketchesData.title', type="textarea" rows=1, placeholder='(Enter a title here)')
  p
    input.note(v-model='sketchesData.note', placeholder='(Add some note here)')
  .constructs
    transition-group(name='constructs-list',
                     enter-active-class='animated flipInX',
                     leave-active-class='animated fadeOut absolute-animation',
                     tag='div')
      .construct(v-for="construct, i in sketchesData.constructs", v-model='sketchesData.constructs[i]',
           @duplicate='duplicateConstruct(i)', @new='newConstruct(i)',
           @delete='deleteConstruct(i)', @moveUp='moveUpConstruct(i)',
           @moveDown='moveDownConstruct(i)',
           is='construct',
           :key='construct.id')
</template>

<script>
import construct from './Construct'
import download from 'downloadjs'
import uuidv1 from 'uuid/v1'
// import json2xlsx from 'json2xlsx-export'
import XLSX from 'xlsx'

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
                id: uuidv1()
              },
              {
                category: 'CDS',
                label: 'gene with a very very long name',
                id: uuidv1()
              },
              {
                category: 'terminator',
                label: 'Term1',
                id: uuidv1()
              }
            ]
          }
        ]
      })
    }
  },
  data () {
    return {
      sketchesData: this.value,
      file: null,
      showFileDialog: false,
      constructIdCounter: 5
    }
  },
  components: {
    construct
  },
  methods: {
    downloadJson () {
      download(JSON.stringify(this.sketchesData, null, ' '), 'design.json')
    },
    downloadXlsx () {
      var data = this.sketchesData
      var workbook = XLSX.utils.book_new()
      var sheet = XLSX.utils.aoa_to_sheet([
        ['field', 'value'],
        ['title', data.title],
        ['note', data.note]
      ])
      XLSX.utils.book_append_sheet(workbook, sheet, 'options')
      var columns = ['label', 'category', 'reversed', 'subscript', 'sublabel', 'bg_color']
      data.constructs.map(function (construct) {
        console.log(construct)
        var sheetData = construct.parts.map(part => (columns.map(function (c) {
          var result = part[c]
          if ((c === 'reversed') && (!result)) {
            console.log('lol')
            return ''
          } else {
            return result
          }
        })))
        sheet = XLSX.utils.aoa_to_sheet([ columns ].concat(sheetData))
        XLSX.utils.book_append_sheet(workbook, sheet, construct.name)
      })
      var wbData = XLSX.write(workbook, { bookType: 'xlsx', bookSST: false, type: 'array' })
      download(wbData, (data.title === '' ? 'design' : data.title) + '.xlsx', 'application/xlsx')
    },
    newConstruct (i) {
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
            id: uuidv1()
          }
        ]
      })
      this.$set(this.sketchesData, 'constructs', newConstructs)
    },
    duplicateConstruct (i) {
      var newConstructs = this.sketchesData.constructs.slice()
      var newConstruct = Object.assign(
        {}, JSON.parse(JSON.stringify(newConstructs[i])),
        {id: uuidv1()}
      )
      newConstructs.splice(i + 1, 0, newConstruct)
      this.$set(this.sketchesData, 'constructs', newConstructs)
    },
    deleteConstruct (i) {
      var newConstructs = this.sketchesData.constructs.slice()
      newConstructs.splice(i, 1)
      this.$set(this.sketchesData, 'constructs', newConstructs)
    },
    moveUpConstruct (i) {
      var indexUp = Math.max(0, i - 1)
      var newConstructs = this.sketchesData.constructs.slice()
      var construct = newConstructs[i]
      newConstructs[i] = newConstructs[indexUp]
      newConstructs[indexUp] = construct
      this.$set(this.sketchesData, 'constructs', newConstructs)
    },
    moveDownConstruct (i) {
      if (i === this.sketchesData.constructs.length - 1) {
        return
      }
      this.moveUpConstruct(i + 1)
    }
  },
  watch: {
    file: {
      deep: true,
      handler (newval) {
        this.showFileDialog = false
        var data = this.sketchesData
        if (newval) {
          var [mimetype, content] = newval.content.split(',')
          if (mimetype.indexOf('json') >= 0) {
            this.sketchesData = JSON.parse(atob(content))
          } else {
            var colorReplacements = {
              blue: '#ECF3FF',
              green: '#DFFFE3',
              red: '#FFEBE9'
            }
            var workbook = XLSX.read(content, {type: 'base64'})
            console.log(workbook)
            workbook.SheetNames.map(function (sheetName) {
              var sheet = workbook.Sheets[sheetName]
              if (sheetName === 'options') {
                XLSX.utils.sheet_to_json(sheet).map(function (o) {
                  if (data[o.field] === '') {
                    data[o.field] = o.value
                  }
                })
              } else {
                var newConstruct = {
                  id: uuidv1(),
                  name: sheetName,
                  note: '',
                  parts: XLSX.utils.sheet_to_json(sheet).map(function (o) {
                    o.bg_color = colorReplacements[o.bg_color] || o.bg_color
                    return Object.assign({id: uuidv1()}, o)
                  })
                }
                data.constructs.push(newConstruct)
              }
            })
          }
          this.file = null
        }
      }
    },
    'sketchesData': {
      deep: true,
      handler (newval) {
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
    text-align: left;
    .file-link {
      display: inline;
      font-size: 1em;
      border: none;
      margin-right: 1em;
    }
    margin-bottom: 3em;

  }
.constructs-list-move {
  transition: transform 1s;
}
// .constructs-list-item {
//   transition: all 1s;
// }
.absolute-animation {
  position: absolute;
  transform: all 0.3s;
}
}
</style>
