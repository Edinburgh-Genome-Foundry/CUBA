<template lang='pug'>
.construct
  hr
  .controls.construct-hover-only
    el-button-group
      el-button(@click="$emit('delete')" size='mini' icon='el-icon-delete' title='Delete')
      el-button(@click="$emit('new')" size='mini' icon='el-icon-plus' title='Insert new')
      el-button(@click="$emit('moveUp')" size='mini' icon='el-icon-arrow-up' title='Move up')
      el-button(@click="$emit('moveDown')" size='mini' icon='el-icon-arrow-down' title='Move down')
      el-button(@click="$emit('duplicate')" size='mini' icon='el-icon-news') Duplicate
  textarea.name(v-model='constructData.name', placeholder='(Name this construct)' rows=1)
  textarea.note(v-model='constructData.note', placeholder='(Add a note)',
                :class="{'construct-hover-only': (constructData.note.length === 0)}")
  .parts
    transition-group.inline-part(name='parts-list',
                     enter-active-class='animated flipInX',
                     leave-active-class='animated zoomOut absolute',
                     tag='div')
      .inline-part(@click="addPart(0)", is='part-adder', key="www")
      .inline-part(v-for='part, i in constructData.parts', :key="part.id")
        part(v-model='constructData.parts[i]' @delete='deletePart(i+1)')
        part-adder(@click="addPart(i+1)")

</template>
<script>

import part from './Part'
import partAdder from './PartAdder'

export default {
  name: 'construct',
  props: {
    value: { default: () => ({
      name: '',
      note: '',
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
    })}
  },
  data () {
    return {
      constructData: this.value,
      counter: 5
    }
  },
  components: {
    part: part,
    'part-adder': partAdder
  },
  methods: {
    addPart (i) {
      var newPartsData = this.constructData.parts.slice()
      this.counter++
      newPartsData.splice(i, 0, {
        category: 'user-defined',
        label: '',
        id: this.counter
      })
      this.$set(this.constructData, 'parts', newPartsData)
    },
    deletePart (i) {
      var newPartsData = this.constructData.parts.slice()
      newPartsData.splice(i - 1, 1)
      this.$set(this.constructData, 'parts', newPartsData)
      console.log('delete', i)
    }
  },
  watch: {
    constructData: {
      deep: true,
      handler (newval) {
        this.$emit('input', newval)
      }
    }
  }
}
</script>
<!-- Add 'scoped' attribute to limit CSS to this component only -->
<style lang='scss' scoped>
.construct {
  hr {
    border: 1px solid #eee;
  }
  margin-top: 4em;
  .name {
    font-size: 1.5em;
    margin-bottom: 0em;
    width: 100%;
  }
  .note {
    width: 100%;
  }
  .controls {
    text-align: right;
    margin-bottom: 1em;
  }
  textarea {
    border: none !important;
    outline: 0 !important;
    border-color: Transparent;
    font-family: 'Raleway';
    font-size: 1em;
    resize: none;

  }
  .construct-hover-only {
    visibility: hidden;
  }
  .sketcher-part.adder /deep/ {
    p.construct-hover-only {
      visibility: hidden;
    }
  }
  &:hover {
    .construct-hover-only {
      visibility: visible;
    }
    .sketcher-part.adder /deep/ {
      p.construct-hover-only {
        visibility: visible;
      }
    }
  }
.inline-part {
  display: inline-block;;
}
.parts-list-move {
  transition: transform 0.5s;
}
.absolute {
  position: absolute;
  transition: all 0.5s;
}
}
</style>
