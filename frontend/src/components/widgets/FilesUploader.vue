https://jsfiddle.net/Linusborg/dzfdctv9/

<template lang="pug">
.dropzone-area(drag-over='handleDragOver', @dragenter='hovering = true',
               @dragleave='hovering = false', :class="{'hovered': hovering}",
               :style="{'height': '100px'}")
    .dropzone-text
      .dropzone-title {{text}}
      .dropzone-info(v-if='help') {{ help }}
    input(type='file', @change='change', multiple='true')
</template>

<script>
export default {
  props: {
    text: {default: 'Drop files here or click to select'},
    help: {default: 'No files too big though :)'},
    filter: {default: () => () => true},
    value: {default: () => ([])}
  },
  data: function () {
    return {
      hovering: false,
      valueMirror: this.value
    }
  },
  methods: {
    change: function (evt) {
      var files = evt.target.files
      var self = this
      self.valueMirror = []
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
            }
          }
          reader.readAsDataURL(files[i])
        }
      }
    }
  },
  watch: {
    valueMirror: function (val) {
      // this.value = val
      this.$emit('input', val)
    }
  }
}
</script>

<style scoped>
.el-icon-upload {
  font-size: 20px;
}

.dropzone-area {
    height: 200px;
    position: relative;
    border: 2px dashed #CBCBCB;
    &.hovered {
        border: 2px dashed #2E94C4;
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
</style>
