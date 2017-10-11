<template lang='pug'>
.construct.animated.flipInX
  hr
  textarea.name(v-model='constructData.name', placeholder='(Name this construct)' rows=1)
  textarea.note(v-model='constructData.note', placeholder='(Add a note)',
                :class="{'construct-hover-only': (constructData.note.length === 0)}")

  .controls.construct-hover-only
    .control(@click="$emit('delete')") <icon name='trash-o'></icon> Delete
    .control(@click="$emit('new')") <icon name='plus-circle'></icon> Insert new
    .control(@click="$emit('duplicate')") <icon name='files-o'></icon> Duplicate
    .control(@click="$emit('moveUp')") <icon name='arrow-up'></icon> Move up
    .control(@click="$emit('moveDown')") <icon name='arrow-down'></icon> Move down
  .parts
    part-adder(@click="addPart(0)")
    span(v-for='part, i in constructData.parts', :key="part.id")
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
  data: function () {
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
    addPart: function (i) {
      var newPartsData = this.constructData.parts.slice()
      this.counter++
      newPartsData.splice(i, 0, {
        category: 'user-defined',
        label: '',
        id: this.counter
      })
      this.$set(this.constructData, 'parts', newPartsData)
      console.log('add', i)
    },
    deletePart: function (i) {
      var newPartsData = this.constructData.parts.slice()
      newPartsData.splice(i - 1, 1)
      this.$set(this.constructData, 'parts', newPartsData)
      console.log('delete', i)
    }
  },
  watch: {
    constructData: {
      deep: true,
      handler: function (newval) {
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
    margin-top: -1em;
    margin-bottom: 1em;


    .control {
      display: inline-block;
      margin-right: 1.5em;
      // width: 6.5em;
      color: #91a0cc;
      cursor: pointer;
      &:hover {
        // font-weight: bold;
        color: #000;
      }
    }
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
}
</style>
