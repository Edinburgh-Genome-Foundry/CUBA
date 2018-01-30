# File uploader

```
usage:
file-uploader(v-model='uploadedFiles')
```

When the user drops/select files, the variable v-model will contain a list
of the form

```
[
  {
    name: 'file1.png',
    content:'data:png64,S0m3GarBaG3dAta...'
  },
  {
    name: 'file2.csv',
    content:'A1, bla, bli\nA2, ble, blu\n...'
  },
  etc.
]
```

From a recipe from here:
https://jsfiddle.net/Linusborg/dzfdctv9/


<template lang="pug">
.files-uploader
  .dropzone-area(@dragenter='hovering = true',
                 @dragover='hovering = true',
                 @dragleave='hovering = false', :class="{'hovered': hovering}",
                 :style="{'height': '100px'}")
      .dropzone-text
        .dropzone-title {{text}}
        .dropzone-info(v-if='help') {{ help }}
      input(type='file', @change='change', :multiple='multiple')
  p(v-if='loadingInProgress && (loaded < toLoad)') Loading file {{loaded + 1}} / {{toLoad}}
  p.selected(v-if='!loadingInProgress && value && ((value.name && !(multiple)) || (value.length && !(value.length === 0)))')

    span Selected:
    span(v-if='value.name') <b> {{value.name}} </b>
    span.selected-file(v-for='file in value' v-else) <b>{{file.name}}</b>
</template>

<script>
export default {
  props: {
    text: {default: 'Drop files here or click to select'},
    help: {default: 'No files too big though :)'},
    filter: {default: () => () => true},
    value: {default: () => ([])},
    displaySelected: {default: true},
    multiple: {default: true}
  },
  data () {
    return {
      hovering: false,
      valueMirror: this.value,
      loadingInProgress: false,
      toLoad: 0,
      loaded: 0
    }
  },
  methods: {
    change (evt) {
      this.hovering = false
      var files = evt.target.files
      var self = this
      self.valueMirror = []

      self.toLoad = files.length
      self.loaded = 0
      self.loadingInProgress = true
      for (var i = 0; i < files.length; i++) {
        var file = files[i]
        let name = file.name
        if (this.filter(files[i])) {
          var reader = new window.FileReader()
          reader.onloadend = function (ev) {
            if (ev.target.readyState === window.FileReader.DONE) {
              self.valueMirror.push({
                name: name,
                content: ev.target.result
              })
              self.loaded++
              if (self.loaded === self.toLoad) {
                self.loadingInProgress = false
              }
            }
          }
          reader.readAsDataURL(files[i])
        }
      }
    }
  },
  watch: {
    valueMirror: {
      deep: true,
      handler (val) {
        this.$emit('input', this.multiple ? val : val[0])
      }
    }
  }
}
</script>

<style lang='scss' scoped>
.el-icon-upload {
  font-size: 20px;
}

.dropzone-area {
    height: 200px;
    margin-bottom:20px;
    position: relative;
    border: 2px dashed #CBCBCB;
    &.hovered {
        border: 2px dashed #2E94C4;
        background-color: rgba(#19749f, 0.03);
        .dropzone-title {
          color: #1975A0;
        }
    }
}

.dropzone-area input {
    position: absolute;
    cursor: pointer;
    top: 0px;
    right: 0;
    bottom: 0;
    left: 0;
    width: 100%;
    height: 100%;
    opacity: 0;
}

.dropzone-text {
    position: absolute;
    top: 50%;
    text-align: center;
    transform: translate(0, -50%);
    width: 100%;
    span {
        display: block;
        font-family: Arial, Helvetica;
        line-height: 1.9;
    }
}

.dropzone-title {
    font-size: 13px;
    color: #787878;
    letter-spacing: 0.4px;
}
.dropzone-info {
    font-size: 13px;
    color: #A8A8A8;
    letter-spacing: 0.4px;
}

.dropzone-button {
    position: absolute;
    top: 10px;
    right: 10px;
    display: none;
}

.dropzone-preview {
    width: 80%;
    position: relative;
    &:hover .dropzone-button {
        display: block;
    }
    img {
      display: block;
      height: auto;
      max-width: 100%;
    }

}

p.selected {
  font-family: 'Inconsolata', Courier;
  font-size: 14px;
  margin-top: -10px;
  .selected-file {
    margin-right: 1em;
  }
}
span.selected-file {
  display: inline-block;
}
</style>
